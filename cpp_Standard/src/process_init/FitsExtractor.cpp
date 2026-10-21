#include "process_init/FitsExtractor.hpp"
#include "InitConfig.hpp"

#include <fitsio.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fqinit {
namespace {
namespace fs = std::filesystem;

constexpr std::size_t kCopyBufferBytes = 16U * 1024U * 1024U;
constexpr std::uint64_t kFitsBlockSize = 2880;
constexpr std::size_t kHeaderMaxBlocks = 8192;  // Bound untrusted header scanning.

struct PlannedImage {
    int hdu_number = 0;
    int expected_ccdnum = 0;
    bool check_ccdnum = false;
    fs::path final_path;
    fs::path staged_path;
    std::uint64_t physical_offset = 0;
    std::uint64_t physical_size = 0;
    bool raw_slice = false;
    bool repair_header = false;
    std::uint64_t header_bytes = 0;
};

struct RawHdu {
    int number = 0;
    std::uint64_t offset = 0;
    std::uint64_t bytes = 0;
    bool is_image = false;
    bool compressed = false;
    int naxis = 0;
    int ccdnum = 0;
    bool has_ccdnum = false;
    bool bad_compression_header = false;
    bool repair_header = false;
    std::uint64_t header_bytes = 0;
};

std::string fitsDiagnostic(int status) {
    char text[FLEN_STATUS] = {};
    fits_get_errstatus(status, text);
    std::ostringstream out;
    out << text;
    char detail[FLEN_ERRMSG] = {};
    while (fits_read_errmsg(detail) != 0) out << "; " << detail;
    return out.str();
}

void closeFits(fitsfile*& file) {
    if (file != nullptr) {
        int status = 0;
        fits_close_file(file, &status);
        file = nullptr;
    }
}

void replaceAll(std::string& text, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    std::size_t pos = 0;
    while ((pos = text.find(from, pos)) != std::string::npos) {
        text.replace(pos, from.size(), to);
        pos += to.size();
    }
}

std::pair<int, std::size_t> imageDatatype(int bitpix) {
    switch (bitpix) {
        case BYTE_IMG:      return {TBYTE, 1};
        case SBYTE_IMG:     return {TSBYTE, 1};
        case SHORT_IMG:     return {TSHORT, 2};
        case USHORT_IMG:    return {TUSHORT, 2};
        case LONG_IMG:      return {TINT, 4};
        case ULONG_IMG:     return {TUINT, 4};
        case LONGLONG_IMG:  return {TLONGLONG, 8};
        case ULONGLONG_IMG: return {TULONGLONG, 8};
        case FLOAT_IMG:     return {TFLOAT, 4};
        case DOUBLE_IMG:    return {TDOUBLE, 8};
        default: throw std::runtime_error("unsupported FITS BITPIX: " + std::to_string(bitpix));
    }
}

std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(' ');
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(' ');
    return s.substr(first, last - first + 1);
}

// These cards are necessary for logical pixels or are read by the main stage.
// The fallback checks TYPE, not arbitrary numerical ranges. A structural
// dimension must nevertheless be positive to allocate/read an image.
bool criticalKeyword(const std::string& name) {
    static const std::set<std::string> keys = {
        "CCDNUM", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
        "BSCALE", "BZERO", "BLANK", "CRPIX1", "CRPIX2",
        "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"
    };
    return keys.count(name) != 0;
}

// Forward declaration: raw FITS numeric syntax checker, implemented below.
bool validFitsNumber(std::string token);

bool nativeNumericHeader(fitsfile* input, const char* keyword, bool integer,
                         bool missing_ok = false) {
    char card[FLEN_CARD] = {};
    int status = 0;
    fits_read_card(input, keyword, card, &status);
    if (status == KEY_NO_EXIST && missing_ok) {
        fits_clear_errmsg();
        return true;
    }
    if (status != 0) {
        const std::string error = fitsDiagnostic(status);
        throw std::runtime_error(std::string("critical FITS header ") + keyword +
                                 " missing or unreadable: " + error);
    }
    std::string c(card, std::min<std::size_t>(std::strlen(card), 80));
    c.resize(80, ' ');
    if (c[8] != '=' || c[9] != ' ')
        return false;
    std::string token = trim(c.substr(10));
    token = trim(token.substr(0, token.find('/')));
    if (token.empty() || token.front() == '\'') return false;
    if (!integer) return validFitsNumber(token);
    std::size_t i = 0;
    if (token[i] == '+' || token[i] == '-') ++i;
    if (i == token.size()) return false;
    for (; i < token.size(); ++i)
        if (token[i] < '0' || token[i] > '9') return false;
    return true;
}

// Shared raw-card grammar helpers, implemented below.
bool printableFitsCard(const std::string& card);
bool invalidCardSyntax(const std::string& card, const std::string& key);
bool structuralOrDecodingKeyword(const std::string& key);

// Only definitively malformed optional cards are removed. Valid unknown
// keywords, including valid custom keywords, are preserved unchanged.
// This runs only on the fallback path; source files are always read-only.
void copyMetadata(fitsfile* input, fitsfile* output, const fs::path& source,
                  int hdu_number) {
    int status = 0;
    int count = 0;
    fits_get_hdrspace(input, &count, nullptr, &status);
    if (status != 0)
        throw std::runtime_error("cannot enumerate source header cards: " + fitsDiagnostic(status));

    for (int i = 1; i <= count; ++i) {
        char card[FLEN_CARD] = {};
        int read_status = 0;
        fits_read_record(input, i, card, &read_status);
        if (read_status != 0) {
            std::cerr << "[init] " << source << " HDU=" << hdu_number
                      << " DROP optional header record " << i << ": "
                      << fitsDiagnostic(read_status) << '\n';
            continue;
        }
        // CFITSIO returns a NUL-terminated card, not necessarily 80 printable
        // bytes in the caller buffer. Restore FITS space-padding before check.
        std::string raw(card, std::min<std::size_t>(std::strlen(card), 80));
        raw.resize(80, ' ');
        const std::string key = trim(raw.substr(0, 8));
        if ((!printableFitsCard(raw) || invalidCardSyntax(raw, key))) {
            if (structuralOrDecodingKeyword(key) || i == 1)
                throw std::runtime_error("malformed critical FITS card " + key);
            std::cerr << "[init] " << source << " HDU=" << hdu_number
                      << " DROP malformed optional header record " << i
                      << " key=" << (key.empty() ? "<blank>" : key) << '\n';
            continue;
        }
        const int cls = fits_get_keyclass(card);
        if (cls <= TYP_CMPRS_KEY || cls == TYP_CKSUM_KEY) continue;
        int write_status = 0;
        fits_write_record(output, card, &write_status);
        if (write_status == 0) continue;
        const std::string why = fitsDiagnostic(write_status);
        if (structuralOrDecodingKeyword(key))
            throw std::runtime_error("cannot copy critical FITS card " + key + ": " + why);
        std::cerr << "[init] " << source << " HDU=" << hdu_number
                  << " DROP CFITSIO-rejected optional header " << key
                  << ": " << why << '\n';
    }
}

// Only required types are verified. No range/finite checks are imposed on
// WCS, pixel scaling, or CCDNUM. NAXIS==2 is enforced at image selection.
// Both Science and DQ require CCDNUM, including when PCA is disabled.
void validateCriticalHeaders(fitsfile* input, ProductKind kind) {
    (void)kind;
    int status = 0;
    int ndim = 0;
    fits_get_img_dim(input, &ndim, &status);
    if (status != 0 || ndim != 2)
        throw std::runtime_error("critical NAXIS is not a two-dimensional image" +
                                 (status ? std::string(": ") + fitsDiagnostic(status) : ""));

    const char* ccd_key = InitConfig::CCDNUM_KEYWORD;
    if (!nativeNumericHeader(input, ccd_key, true))
        throw std::runtime_error("critical CCDNUM must be stored as a FITS integer, not a string/float");
    int ccdnum = 0;
    status = 0;
    fits_read_key(input, TINT, ccd_key, &ccdnum, nullptr, &status);
    if (status != 0)
        throw std::runtime_error(std::string("critical CCDNUM must be integer: ") +
                                 fitsDiagnostic(status));

    for (const char* key : {"CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                            "CD1_1", "CD1_2", "CD2_1", "CD2_2"}) {
        if (kind != ProductKind::Science) break;
        if (!nativeNumericHeader(input, key, false))
            throw std::runtime_error(std::string("critical science WCS ") + key +
                                     " must be stored as a FITS number");
        double value = 0.0;
        status = 0;
        fits_read_key(input, TDOUBLE, key, &value, nullptr, &status);
        if (status != 0)
            throw std::runtime_error(std::string("critical science WCS ") + key +
                                     " must be numeric: " + fitsDiagnostic(status));
    }
    for (const char* key : {"BSCALE", "BZERO"}) {
        if (!nativeNumericHeader(input, key, false, true))
            throw std::runtime_error(std::string("pixel scaling ") + key +
                                     " must be stored as a FITS number");
        double value = 0.0;
        status = 0;
        fits_read_key(input, TDOUBLE, key, &value, nullptr, &status);
        if (status == KEY_NO_EXIST) { fits_clear_errmsg(); continue; }
        if (status != 0)
            throw std::runtime_error(std::string("pixel scaling ") + key +
                                     " must be numeric: " + fitsDiagnostic(status));
    }
    if (!nativeNumericHeader(input, "BLANK", true, true))
        throw std::runtime_error("BLANK must be stored as a FITS integer");
    LONGLONG blank = 0;
    status = 0;
    fits_read_key(input, TLONGLONG, "BLANK", &blank, nullptr, &status);
    if (status != 0 && status != KEY_NO_EXIST)
        throw std::runtime_error("BLANK must be integer: " + fitsDiagnostic(status));
    if (status == KEY_NO_EXIST) fits_clear_errmsg();
}

// Existing production extraction is unchanged for well-formed FITS files.
// Only a failure from copying header records triggers the relaxed path.
struct HeaderCopyFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

void copyCurrentImageOriginal(fitsfile* input, const fs::path& output_path) {
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_get_img_param(input, 9, &bitpix, &naxis, naxes, &status);
    if (status != 0)
        throw std::runtime_error("cannot read logical image parameters: " + fitsDiagnostic(status));
    if (naxis != 2)
        throw std::runtime_error("initializer only extracts two-dimensional image HDUs");

    fitsfile* output = nullptr;
    try {
        const std::string create_name = "!" + output_path.string();
        fits_create_file(&output, create_name.c_str(), &status);
        if (status == 0) fits_create_img(output, bitpix, naxis, naxes, &status);
        if (status != 0)
            throw std::runtime_error("cannot create uncompressed output: " + fitsDiagnostic(status));
        int key_count = 0;
        fits_get_hdrspace(input, &key_count, nullptr, &status);
        for (int key_index = 1; key_index <= key_count && status == 0; ++key_index) {
            char card[FLEN_CARD] = {};
            fits_read_record(input, key_index, card, &status);
            if (status != 0) break;
            const int key_class = fits_get_keyclass(card);
            if (key_class > TYP_CMPRS_KEY && key_class != TYP_CKSUM_KEY)
                fits_write_record(output, card, &status);
        }
        if (status != 0)
            throw HeaderCopyFailure("cannot copy FITS header cards: " + fitsDiagnostic(status));

        const auto [datatype, bytes_per_pixel] = imageDatatype(bitpix);
        LONGLONG total_pixels = 1;
        for (int axis = 0; axis < naxis; ++axis) {
            if (naxes[axis] <= 0)
                throw std::runtime_error("image has a non-positive axis length");
            if (total_pixels > LLONG_MAX / naxes[axis])
                throw std::runtime_error("image pixel count overflow");
            total_pixels *= static_cast<LONGLONG>(naxes[axis]);
        }
        const LONGLONG chunk_pixels = std::max<LONGLONG>(
            1, static_cast<LONGLONG>(kCopyBufferBytes / bytes_per_pixel));
        const std::size_t buffer_bytes = static_cast<std::size_t>(chunk_pixels) * bytes_per_pixel;
        std::vector<double> aligned_buffer((buffer_bytes + sizeof(double) - 1U) / sizeof(double));
        double null_value = 0.0;
        int any_null = 0;
        fits_set_bscale(input, 1.0, 0.0, &status);
        fits_set_bscale(output, 1.0, 0.0, &status);
        LONGLONG first_pixel = 1;
        LONGLONG remaining = total_pixels;
        while (remaining > 0 && status == 0) {
            const LONGLONG count = std::min(remaining, chunk_pixels);
            fits_read_img(input, datatype, first_pixel, count, &null_value,
                          aligned_buffer.data(), &any_null, &status);
            if (status == 0)
                fits_write_img(output, datatype, first_pixel, count, aligned_buffer.data(), &status);
            first_pixel += count;
            remaining -= count;
        }
        int close_status = 0;
        fits_close_file(output, &close_status);
        output = nullptr;
        if (status == 0 && close_status != 0) status = close_status;
        if (status != 0)
            throw std::runtime_error("cannot stream FITS pixels: " + fitsDiagnostic(status));
    } catch (...) {
        closeFits(output);
        throw;
    }
}

void copyCurrentImageRelaxed(fitsfile* input, const fs::path& output_path,
                             const fs::path& source, int hdu_number, ProductKind kind) {
    int status = 0;
    int bitpix = 0;
    int naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_get_img_param(input, 9, &bitpix, &naxis, naxes, &status);
    if (status != 0) throw std::runtime_error("cannot read image parameters: " + fitsDiagnostic(status));
    if (naxis != 2 || naxes[0] <= 0 || naxes[1] <= 0)
        throw std::runtime_error("invalid non-positive or non-two-dimensional logical image");
    const auto [datatype, bytes_per_pixel] = imageDatatype(bitpix);
    if (static_cast<unsigned long long>(naxes[0]) >
        static_cast<unsigned long long>(LLONG_MAX) / static_cast<unsigned long long>(naxes[1]))
        throw std::runtime_error("pixel count overflow");
    const LONGLONG total = static_cast<LONGLONG>(naxes[0]) * naxes[1];
    validateCriticalHeaders(input, kind);

    fitsfile* output = nullptr;
    try {
        const std::string filename = "!" + output_path.string();
        fits_create_file(&output, filename.c_str(), &status);
        if (status != 0) throw std::runtime_error("cannot create FITS: " + fitsDiagnostic(status));
        fits_create_img(output, bitpix, 2, naxes, &status);
        if (status != 0) throw std::runtime_error("cannot create FITS image: " + fitsDiagnostic(status));
        copyMetadata(input, output, source, hdu_number);
        const LONGLONG chunk = std::max<LONGLONG>(1,
            static_cast<LONGLONG>(kCopyBufferBytes / bytes_per_pixel));
        std::vector<double> buffer((static_cast<std::size_t>(chunk) *
                                    bytes_per_pixel + sizeof(double) - 1) / sizeof(double));
        fits_set_bscale(input, 1.0, 0.0, &status);
        if (status == 0) fits_set_bscale(output, 1.0, 0.0, &status);
        if (status != 0) throw std::runtime_error("cannot set FITS raw pixel scale: " + fitsDiagnostic(status));
        double null_value = 0.0;
        int any_null = 0;
        for (LONGLONG first = 1; first <= total;) {
            const LONGLONG n = std::min(chunk, total - first + 1);
            fits_read_img(input, datatype, first, n, &null_value,
                          buffer.data(), &any_null, &status);
            if (status != 0) throw std::runtime_error("FITS pixel read failed: " + fitsDiagnostic(status));
            fits_write_img(output, datatype, first, n, buffer.data(), &status);
            if (status != 0) throw std::runtime_error("FITS pixel write failed: " + fitsDiagnostic(status));
            first += n;
        }
        fits_close_file(output, &status);
        output = nullptr;
        if (status != 0) throw std::runtime_error("FITS output close failed: " + fitsDiagnostic(status));
    } catch (...) {
        closeFits(output);
        throw;
    }
}

bool validateOutput(const PlannedImage& plan, std::string& error) {
    fitsfile* file = nullptr;
    int status = 0;
    fits_open_file(&file, plan.final_path.string().c_str(), READONLY, &status);
    if (status != 0) {
        error = "cannot open existing output " + plan.final_path.string() + ": " + fitsDiagnostic(status);
        closeFits(file);
        return false;
    }
    int naxis = 0;
    fits_get_img_dim(file, &naxis, &status);
    const int compressed = fits_is_compressed_image(file, &status);
    if (status == 0 && (naxis != 2 || compressed != 0)) {
        error = "existing output is not an uncompressed two-dimensional FITS image: "
                + plan.final_path.string();
        closeFits(file);
        return false;
    }
    if (status == 0 && plan.check_ccdnum) {
        try {
            if (!nativeNumericHeader(file, InitConfig::CCDNUM_KEYWORD, true))
                status = BAD_KEYCHAR;
        } catch (const std::exception& e) {
            error = std::string("invalid existing CCDNUM: ") + e.what();
            closeFits(file);
            return false;
        }
        int ccdnum = 0;
        if (status == 0)
            fits_read_key(file, TINT, InitConfig::CCDNUM_KEYWORD, &ccdnum, nullptr, &status);
        if (status == 0 && ccdnum != plan.expected_ccdnum) {
            error = "existing output has the wrong " + std::string(InitConfig::CCDNUM_KEYWORD)
                    + ": " + plan.final_path.string();
            closeFits(file);
            return false;
        }
    }
    if (status != 0) {
        error = "cannot validate existing output " + plan.final_path.string() + ": " + fitsDiagnostic(status);
        closeFits(file);
        return false;
    }
    closeFits(file);
    return true;
}

std::vector<PlannedImage> planImages(fitsfile* input, const fs::path& source,
                                     ProductKind kind, const fs::path& final_dir,
                                     const fs::path& staging_dir, int& skipped) {
    int status = 0, nhdu = 0;
    fits_get_num_hdus(input, &nhdu, &status);
    if (status != 0) throw std::runtime_error("cannot count FITS HDUs: " + fitsDiagnostic(status));
    const std::string stem = kind == ProductKind::Science ? archiveStem(source) : dqOutputStem(source);
    int science_index = 0;
    std::set<std::string> unique_names;
    std::vector<PlannedImage> plans;
    for (int n = 1; n <= nhdu; ++n) {
        int type = 0, s = 0;
        fits_movabs_hdu(input, n, &type, &s);
        if (s != 0) {
            std::cerr << "[init] " << source << " HDU=" << n
                      << " header navigation failed: " << fitsDiagnostic(s) << '\n';
            continue;
        }
        if (type != IMAGE_HDU) continue;
        int naxis = 0;
        s = 0;
        fits_get_img_dim(input, &naxis, &s);
        if (s != 0) {
            std::cerr << "[init] " << source << " HDU=" << n
                      << " dimensions failed: " << fitsDiagnostic(s) << '\n';
            continue;
        }
        if (naxis != 2) continue;
        // Reserve the Science ordinal before checking CCDNUM. This keeps
        // original filenames stable if an earlier image HDU is skipped.
        const int science_number = kind == ProductKind::Science ? ++science_index : 0;
        int ccdnum = 0;
        s = 0;
        bool native_ccd = false;
        try {
            native_ccd = nativeNumericHeader(input, InitConfig::CCDNUM_KEYWORD, true);
        } catch (const std::exception& error) {
            std::cerr << "[init] " << source << " HDU=" << n
                      << " invalid CCDNUM: " << error.what() << '\n';
        }
        if (!native_ccd) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << n
                      << " SKIP: CCDNUM is not stored as a FITS integer\n";
            fits_clear_errmsg();
            continue;
        }
        fits_read_key(input, TINT, InitConfig::CCDNUM_KEYWORD, &ccdnum, nullptr, &s);
        if (s != 0) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << n
                      << " SKIP: CCDNUM cannot be converted to int: " << fitsDiagnostic(s) << '\n';
            continue;
        }
        const int number = kind == ProductKind::DqMask ? ccdnum : science_number;
        const bool check = true;
        const std::string filename = stem + "_" + std::to_string(number) + ".fits";
        if (!unique_names.insert(filename).second) {
            std::cerr << "[init] " << source << " HDU=" << n
                      << " skipped: duplicate output " << filename << '\n';
            continue;
        }
        plans.push_back({n, ccdnum, check, final_dir / filename, staging_dir / filename});
    }
    if (plans.empty()) throw std::runtime_error("archive contains no extractable two-dimensional image HDUs");
    return plans;
}

// Read raw FITS headers without invoking CFITSIO. Used to find surviving HDUs
// if CFITSIO rejects another extension during fits_get_num_hdus(). Each valid
// extension is extracted from a small temporary FITS containing only that HDU.
// This avoids a bad extension's *logical* header poisoning every subsequent HDU.
// Physical FITS structure must still be parseable to find later HDU offsets.

// Only header cards are inspected; pixel data is never pre-read. The cards
// are 80 bytes and contain printable 7-bit FITS ASCII. Unsafe optional cards
// can be rewritten *in a temporary copy*, without mutating the source archive.
bool validFitsKeyword(const std::string& key) {
    if (key.empty() || key.size() > 8) return false;
    for (unsigned char c : key)
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-'))
            return false;
    return true;
}

bool structuralOrDecodingKeyword(const std::string& key) {
    if (criticalKeyword(key)) return true;
    if (key == "SIMPLE" || key == "XTENSION" || key == "PCOUNT" ||
        key == "GCOUNT" || key == "TFIELDS" || key == "THEAP" ||
        key == "ZIMAGE" || key == "ZCMPTYPE" || key == "ZBITPIX" ||
        key == "ZNAXIS" || key == "ZQUANTIZ" || key == "ZDITHER0") return true;
    if (key.rfind("ZNAXIS", 0) == 0 || key.rfind("ZTILE", 0) == 0 ||
        key.rfind("ZNAME", 0) == 0 || key.rfind("ZVAL", 0) == 0 ||
        key.rfind("TTYPE", 0) == 0 || key.rfind("TFORM", 0) == 0 ||
        key.rfind("TSCAL", 0) == 0 || key.rfind("TZERO", 0) == 0 ||
        key.rfind("TNULL", 0) == 0) return true;
    return false;
}

bool printableFitsCard(const std::string& card) {
    return std::all_of(card.begin(), card.end(), [](unsigned char c) {
        return c >= 32 && c <= 126;
    });
}

// FITS numeric syntax (not a semantic range check). Avoid strtod here:
// it accepts non-FITS tokens such as nan, inf and hexadecimal floats.
bool validFitsNumber(std::string token) {
    token = trim(token);
    if (token.empty()) return false;
    std::size_t i = 0;
    if (token[i] == '+' || token[i] == '-') ++i;
    std::size_t digits = 0;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
        ++i; ++digits;
    }
    if (i < token.size() && token[i] == '.') {
        ++i;
        while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
            ++i; ++digits;
        }
    }
    if (digits == 0) return false;
    if (i < token.size() && (token[i] == 'E' || token[i] == 'D' ||
                              token[i] == 'e' || token[i] == 'd')) {
        ++i;
        if (i < token.size() && (token[i] == '+' || token[i] == '-')) ++i;
        const std::size_t exp_begin = i;
        while (i < token.size() && token[i] >= '0' && token[i] <= '9') ++i;
        if (exp_begin == i) return false;
    }
    return i == token.size();
}

// Deliberately conservative FITS card syntax check. It does NOT inspect
// scientific meaning or an optional keyword's expected numeric type.
// HIERARCH/CONTINUE and commentary cards are preserved when printable.
bool invalidCardSyntax(const std::string& card, const std::string& key) {
    if (card.size() != 80) return true;
    if (key.empty()) return false;  // FITS blank-keyword commentary.
    if (!validFitsKeyword(key)) return true;
    if (key == "COMMENT" || key == "HISTORY" || key == "CONTINUE" ||
        key == "HIERARCH" || key == "END") return false;
    if (card[8] != '=' || card[9] != ' ') return true;
    const std::string value = trim(card.substr(10));
    if (value.empty()) return true;
    if (value.front() == '\'') {
        // Quoted FITS strings escape a quote by doubling it. A trailing
        // comment, if present, starts only after the closing quote.
        bool closed = false;
        for (std::size_t j = 11; j < 80; ++j) {
            if (card[j] != '\'') continue;
            if (j + 1 < 80 && card[j + 1] == '\'') { ++j; continue; }
            closed = true;
            break;
        }
        return !closed;
    }
    std::string token = trim(value.substr(0, value.find('/')));
    if (token == "T" || token == "F") return false;
    // FITS complex literals use two comma-separated real numbers.
    if (token.size() >= 5 && token.front() == '(' && token.back() == ')') {
        std::size_t comma = token.find(',');
        if (comma == std::string::npos || comma == 1 || comma + 1 >= token.size() - 1)
            return true;
        return !validFitsNumber(token.substr(1, comma - 1)) ||
               !validFitsNumber(token.substr(comma + 1, token.size() - comma - 2));
    }
    return !validFitsNumber(token);
}

bool needsRawRemoval(const std::string& card, const std::string& key) {
    return !printableFitsCard(card) || invalidCardSyntax(card, key);
}

std::map<std::string, std::string> readRawHeader(std::ifstream& in,
                                                 std::uint64_t start,
                                                 std::uint64_t& header_bytes,
                                                 bool& needs_repair, bool& bad_critical, bool primary_hdu) {
    std::map<std::string, std::string> hdr;
    in.clear();
    in.seekg(static_cast<std::streamoff>(start), std::ios::beg);
    if (!in) throw std::runtime_error("cannot seek FITS header");
    bool end = false;
    std::size_t card_index = 0;
    for (std::size_t b = 0; b < kHeaderMaxBlocks && !end; ++b) {
        std::array<char, 2880> block{};
        in.read(block.data(), block.size());
        if (!in) throw std::runtime_error("truncated FITS header");
        for (std::size_t i = 0; i < 36; ++i) {
            const std::string card(block.data() + i * 80, 80);
            const std::string key = trim(card.substr(0, 8));
            if (card_index == 0 && key != (primary_hdu ? "SIMPLE" : "XTENSION"))
                bad_critical = true;
            ++card_index;
            if (key == "END") {
                if (card != std::string("END") + std::string(77, ' '))
                    needs_repair = true;
                end = true;
                break;
            }
            if (needsRawRemoval(card, key)) {
                if (structuralOrDecodingKeyword(key)) bad_critical = true;
                else needs_repair = true;
            }
            if (card[8] != '=' || key.empty()) continue;
            if (hdr.count(key) && structuralOrDecodingKeyword(key))
                bad_critical = true;
            std::string value = trim(card.substr(10));
            if (!value.empty() && value.front() == '\'') {
                const std::size_t q = value.find('\'', 1);
                if (q != std::string::npos) value = value.substr(1, q - 1);
            } else {
                value = trim(value.substr(0, value.find('/')));
            }
            if (!hdr.count(key)) hdr[key] = value;
        }
        if (end) { header_bytes = (b + 1) * kFitsBlockSize; return hdr; }
    }
    throw std::runtime_error("FITS END card not found within limit");
}

bool rawInteger(const std::map<std::string, std::string>& h,
                const std::string& key, long long& out) {
    auto it = h.find(key);
    if (it == h.end()) return false;
    try {
        std::size_t p = 0;
        long long value = std::stoll(it->second, &p);
        if (p != it->second.size()) return false;
        out = value;
        return true;
    } catch (...) { return false; }
}

bool mulChecked(std::uint64_t x, std::uint64_t y, std::uint64_t& result) {
    if (y && x > UINT64_MAX / y) return false;
    result = x * y;
    return true;
}

std::vector<RawHdu> scanRawHdus(const fs::path& source, bool& incomplete) {
    const auto file_size = fs::file_size(source);
    std::ifstream in(source, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open raw FITS file");
    std::vector<RawHdu> hdus;
    std::uint64_t pos = 0;
    incomplete = false;
    while (pos < file_size) {
      try {
        std::uint64_t header_bytes = 0;
        bool needs_repair = false, bad_critical = false;
        const auto h = readRawHeader(in, pos, header_bytes, needs_repair, bad_critical, pos == 0);
        long long bitpix = 0, naxis = 0, pcount = 0, gcount = 1;
        if (!rawInteger(h, "BITPIX", bitpix) || !rawInteger(h, "NAXIS", naxis) ||
            naxis < 0 || naxis > 9 ||
            (bitpix != 8 && bitpix != 16 && bitpix != 32 && bitpix != 64 &&
             bitpix != -32 && bitpix != -64))
            throw std::runtime_error("unreadable physical BITPIX/NAXIS at byte " + std::to_string(pos));
        (void)rawInteger(h, "PCOUNT", pcount);
        (void)rawInteger(h, "GCOUNT", gcount);
        if (pcount < 0 || gcount <= 0) throw std::runtime_error("invalid physical PCOUNT/GCOUNT");
        std::uint64_t elements = naxis == 0 ? 0 : 1;
        for (int a = 1; a <= naxis; ++a) {
            long long size = 0;
            if (!rawInteger(h, "NAXIS" + std::to_string(a), size) || size < 0 ||
                !mulChecked(elements, static_cast<std::uint64_t>(size), elements))
                throw std::runtime_error("invalid physical axis at byte " + std::to_string(pos));
        }
        std::uint64_t data = 0;
        // For tables, PCOUNT is additional heap bytes, not pixels.
        const bool table = h.count("XTENSION") &&
                           (h.at("XTENSION") == "BINTABLE" || h.at("XTENSION") == "TABLE");
        if (table) {
            if (UINT64_MAX - elements < static_cast<std::uint64_t>(pcount) ||
                !mulChecked(elements + static_cast<std::uint64_t>(pcount),
                            static_cast<std::uint64_t>(gcount), data))
                throw std::runtime_error("physical table size overflow");
        } else if (naxis != 0) {
            std::uint64_t sum = elements;
            if (UINT64_MAX - sum < static_cast<std::uint64_t>(pcount) ||
                !mulChecked(sum + static_cast<std::uint64_t>(pcount),
                            static_cast<std::uint64_t>(gcount), sum) ||
                !mulChecked(sum, static_cast<std::uint64_t>(std::abs(bitpix) / 8), data))
                throw std::runtime_error("physical image size overflow");
        }
        if (data > UINT64_MAX - (kFitsBlockSize - 1))
            throw std::runtime_error("FITS size overflow");
        const std::uint64_t padded = (data + kFitsBlockSize - 1) / kFitsBlockSize * kFitsBlockSize;
        if (header_bytes > file_size - pos || padded > file_size - pos - header_bytes)
            throw std::runtime_error("truncated FITS extension at byte " + std::to_string(pos));
        RawHdu record;
        record.number = static_cast<int>(hdus.size()) + 1;
        record.offset = pos;
        record.bytes = header_bytes + padded;
        record.header_bytes = header_bytes;
        record.repair_header = needs_repair;
        record.bad_compression_header = bad_critical;
        record.compressed = table && h.count("ZIMAGE") && h.at("ZIMAGE") == "T";
        if (record.compressed) {
            long long ndim = 0;
            const bool has_dim = rawInteger(h, "ZNAXIS", ndim);
            record.naxis = has_dim && ndim >= 0 && ndim <= 9 ? static_cast<int>(ndim) : 0;
            record.is_image = true;
            if (!has_dim || ndim != 2) record.bad_compression_header = true;
            for (int axis = 1; axis <= record.naxis; ++axis) {
                long long nx = 0, tile = axis == 1 ? 0 : 1;
                const bool got = rawInteger(h, "ZNAXIS" + std::to_string(axis), nx);
                if (!got) { record.bad_compression_header = true; continue; }
                const std::string tilekey = "ZTILE" + std::to_string(axis);
                if (h.count(tilekey)) {
                    if (!rawInteger(h, tilekey, tile)) record.bad_compression_header = true;
                } else if (axis == 1) {
                    tile = nx;
                }
                if (nx <= 0 || tile <= 0) record.bad_compression_header = true;
            }
            if (h.count("ZCMPTYPE") &&
                (h.at("ZCMPTYPE") == "RICE_1" || h.at("ZCMPTYPE") == "RICE_ONE")) {
                long long bs = 0;
                if (!rawInteger(h, "ZVAL1", bs) || bs <= 0)
                    record.bad_compression_header = true;
            }
        } else {
            record.naxis = static_cast<int>(naxis);
            record.is_image = !table;
        }
        long long ccd = 0;
        record.has_ccdnum = rawInteger(h, InitConfig::CCDNUM_KEYWORD, ccd) &&
                            ccd >= INT_MIN && ccd <= INT_MAX;
        if (record.has_ccdnum) record.ccdnum = static_cast<int>(ccd);
        hdus.push_back(record);
        pos += record.bytes;
      } catch (const std::exception& e) {
        if (hdus.empty()) throw;
        incomplete = true;
        std::cerr << "[init] " << source << " physical FITS layout fails after "
                  << hdus.size() << " HDUs: " << e.what()
                  << "; salvaging preceding HDUs only\n";
        break;
      }
    }
    return hdus;
}

std::vector<PlannedImage> rawPlans(const fs::path& source, ProductKind kind,
                                   const fs::path& final_dir, const fs::path& staging_dir,
                                   const std::vector<RawHdu>& hdus, int& invalid) {
    const std::string stem = kind == ProductKind::Science ? archiveStem(source) : dqOutputStem(source);
    int science_index = 0;
    std::set<std::string> names;
    std::vector<PlannedImage> plans;
    for (const auto& hdu : hdus) {
        if (!hdu.is_image) continue;
        // A damaged compressed image may have an unreadable ZNAXIS. Reserve
        // its position rather than silently shifting later CCD filenames.
        if (hdu.naxis != 2 && !(hdu.compressed && hdu.bad_compression_header)) continue;
        const int ordinal = kind == ProductKind::Science ? ++science_index : 0;
        if (hdu.bad_compression_header || hdu.naxis != 2) {
            ++invalid;
            std::cerr << "[init] " << source << " HDU=" << hdu.number
                      << " SKIP: unsafe critical FITS header or compression parameters\n";
            continue;
        }
        const bool check = true;
        const int number = kind == ProductKind::DqMask ? hdu.ccdnum : ordinal;
        if (!hdu.has_ccdnum) {
            ++invalid;
            std::cerr << "[init] " << source << " HDU=" << hdu.number
                      << " SKIP: critical CCDNUM missing or not integer\n";
            continue;
        }
        const std::string filename = stem + "_" + std::to_string(number) + ".fits";
        if (!names.insert(filename).second) {
            ++invalid;
            std::cerr << "[init] " << source << " HDU=" << hdu.number
                      << " SKIP: duplicate CCD filename " << filename << '\n';
            continue;
        }
        plans.push_back({hdu.number, hdu.ccdnum, check,
                         final_dir / filename, staging_dir / filename,
                         hdu.offset, hdu.bytes, true, hdu.repair_header, hdu.header_bytes});
    }
    return plans;
}

void copyRawRange(const fs::path& source, std::ofstream& out,
                  std::uint64_t offset, std::uint64_t count) {
    std::ifstream in(source, std::ios::binary);
    if (!in) throw std::runtime_error("cannot reopen source raw FITS");
    in.seekg(static_cast<std::streamoff>(offset));
    if (!in) throw std::runtime_error("cannot seek source raw FITS");
    std::vector<char> buffer(1024 * 1024);
    while (count) {
        auto n = static_cast<std::streamsize>(std::min<std::uint64_t>(buffer.size(), count));
        in.read(buffer.data(), n);
        if (in.gcount() != n) throw std::runtime_error("raw FITS slice truncated");
        out.write(buffer.data(), n);
        if (!out) throw std::runtime_error("raw FITS slice write failed");
        count -= static_cast<std::uint64_t>(n);
    }
}


// Rebuild only the TEMPORARY isolated input header. Corrupt optional cards
// are REMOVED and the remaining cards are compacted before END. The physical
// image/table data after the original header are copied unchanged.
void copySanitizedRawHdu(const fs::path& source, std::ofstream& out,
                         std::uint64_t offset, std::uint64_t header_bytes,
                         std::uint64_t total_bytes, int hdu_number) {
    if (header_bytes > total_bytes || header_bytes % kFitsBlockSize != 0)
        throw std::runtime_error("invalid temporary FITS header boundaries");
    std::ifstream in(source, std::ios::binary);
    if (!in) throw std::runtime_error("cannot reopen source for temporary header repair");
    in.seekg(static_cast<std::streamoff>(offset));
    if (!in) throw std::runtime_error("cannot seek source HDU header");

    std::vector<std::array<char, 80>> kept;
    bool found_end = false;
    std::size_t card_index = 0;
    for (std::uint64_t done = 0; done < header_bytes && !found_end;
         done += kFitsBlockSize) {
        std::array<char, 2880> block{};
        in.read(block.data(), block.size());
        if (!in) throw std::runtime_error("truncated source FITS header");
        for (std::size_t j = 0; j < 36; ++j, ++card_index) {
            const std::string card(block.data() + j * 80, 80);
            const std::string key = trim(card.substr(0, 8));
            if (card_index == 0 && key != (hdu_number == 1 ? "SIMPLE" : "XTENSION"))
                throw std::runtime_error("damaged first mandatory HDU header card");
            if (key == "END") { found_end = true; break; }
            if (needsRawRemoval(card, key)) {
                if (card_index == 0 || structuralOrDecodingKeyword(key))
                    throw std::runtime_error("unsafe mandatory/decoding card: " + key);
                std::cerr << "[init] " << source << " HDU=" << hdu_number
                          << " DROP malformed optional source card " << card_index + 1
                          << " key=" << (key.empty() ? "<blank>" : key) << '\n';
                continue;
            }
            std::array<char, 80> copy{};
            std::copy(card.begin(), card.end(), copy.begin());
            kept.push_back(copy);
        }
    }
    if (!found_end) throw std::runtime_error("temporary FITS header missing END");

    std::array<char, 80> end{};
    end.fill(' ');
    std::memcpy(end.data(), "END", 3);
    kept.push_back(end);
    const std::size_t blocks = (kept.size() + 35) / 36;
    for (std::size_t bi = 0; bi < blocks; ++bi) {
        std::array<char, 2880> block{};
        block.fill(' ');
        for (std::size_t j = 0; j < 36 && bi * 36 + j < kept.size(); ++j)
            std::copy(kept[bi * 36 + j].begin(), kept[bi * 36 + j].end(),
                      block.data() + j * 80);
        out.write(block.data(), block.size());
        if (!out) throw std::runtime_error("cannot write compacted temporary FITS header");
    }
    copyRawRange(source, out, offset + header_bytes, total_bytes - header_bytes);
}

// A fresh empty primary is required for a physical IMAGE extension/BINTABLE.
void writeEmptyPrimary(std::ofstream& out) {
    std::array<char, 2880> block{};
    block.fill(' ');
    const std::array<const char*, 5> cards = {
        "SIMPLE  =                    T", "BITPIX  =                    8",
        "NAXIS   =                    0", "EXTEND  =                    T", "END"
    };
    for (std::size_t i = 0; i < cards.size(); ++i) {
        std::memcpy(block.data() + i * 80, cards[i], std::strlen(cards[i]));
    }
    out.write(block.data(), block.size());
}

void extractOne(const fs::path& source, const PlannedImage& plan,
                ProductKind kind, const fs::path& staging_dir,
                fitsfile* original) {
    if (!plan.raw_slice) {
        int type = 0, status = 0;
        fits_movabs_hdu(original, plan.hdu_number, &type, &status);
        if (status != 0 || type != IMAGE_HDU)
            throw std::runtime_error("cannot navigate to planned image: " + fitsDiagnostic(status));
        validateCriticalHeaders(original, kind);
        try {
            copyCurrentImageOriginal(original, plan.staged_path);
        } catch (const HeaderCopyFailure& error) {
            std::cerr << "[init] " << source << " HDU=" << plan.hdu_number
                      << " normal header copy failed; trying temporary relaxed header: "
                      << error.what() << '\n';
            fits_clear_errmsg();
            copyCurrentImageRelaxed(original, plan.staged_path, source, plan.hdu_number, kind);
        }
        return;
    }
    const fs::path slice = staging_dir /
                           (".raw_hdu_" + std::to_string(plan.hdu_number) + ".fits");
    fitsfile* input = nullptr;
    try {
        {
            std::ofstream out(slice, std::ios::binary | std::ios::trunc);
            if (!out) throw std::runtime_error("cannot create temporary isolated HDU");
            if (plan.hdu_number != 1) writeEmptyPrimary(out);
            if (plan.repair_header)
                copySanitizedRawHdu(source, out, plan.physical_offset,
                                    plan.header_bytes, plan.physical_size, plan.hdu_number);
            else
                copyRawRange(source, out, plan.physical_offset, plan.physical_size);
            out.close();
            if (!out) throw std::runtime_error("cannot finalize isolated HDU");
        }
        int status = 0, type = 0;
        fits_open_file(&input, slice.string().c_str(), READONLY, &status);
        if (status == 0) fits_movabs_hdu(input, plan.hdu_number == 1 ? 1 : 2, &type, &status);
        if (status != 0 || type != IMAGE_HDU)
            throw std::runtime_error("isolated HDU cannot be opened: " + fitsDiagnostic(status));
        validateCriticalHeaders(input, kind);
        if (plan.repair_header) {
            copyCurrentImageRelaxed(input, plan.staged_path, source, plan.hdu_number, kind);
        } else {
            try {
                copyCurrentImageOriginal(input, plan.staged_path);
            } catch (const HeaderCopyFailure& error) {
                std::cerr << "[init] " << source << " HDU=" << plan.hdu_number
                          << " header copy failed in isolated HDU; retry relaxed: "
                          << error.what() << '\n';
                fits_clear_errmsg();
                copyCurrentImageRelaxed(input, plan.staged_path, source, plan.hdu_number, kind);
            }
        }
        closeFits(input);
        std::error_code ec;
        fs::remove(slice, ec);
    } catch (...) {
        closeFits(input);
        std::error_code ec;
        fs::remove(slice, ec);
        throw;
    }
}

std::vector<PlannedImage> commitImages(const std::vector<PlannedImage>& plans,
                                       const fs::path& source, int& skipped) {
    std::vector<PlannedImage> committed;
    for (const auto& plan : plans) {
        if (std::rename(plan.staged_path.c_str(), plan.final_path.c_str()) != 0) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << plan.hdu_number
                      << " SKIP: cannot commit " << plan.final_path
                      << ": " << std::strerror(errno) << '\n';
            continue;
        }
        committed.push_back(plan);
    }
    return committed;
}

} // namespace

std::string archiveStem(const fs::path& source) {
    const std::string name = source.filename().string();
    const std::size_t suffix = std::strlen(InitConfig::ARCHIVE_SUFFIX);
    if (name.size() <= suffix ||
        name.compare(name.size() - suffix, suffix, InitConfig::ARCHIVE_SUFFIX) != 0)
        throw std::runtime_error("invalid archive suffix: " + source.string());
    return name.substr(0, name.size() - suffix);
}

std::string dqOutputStem(const fs::path& source) {
    std::string stem = archiveStem(source);
    replaceAll(stem, InitConfig::DQ_STEM_REPLACE_FROM, InitConfig::DQ_STEM_REPLACE_TO);
    return stem;
}

ExtractionResult extractArchive(const fs::path& source, ProductKind kind,
                                const fs::path& final_directory,
                                const fs::path& staging_directory,
                                ExistingPolicy policy) {
    ExtractionResult result;
    fitsfile* input = nullptr;
    try {
        // The raw scanner also checks compressed header dimensions before any
        // CFITSIO HDU scan. For suspicious headers, never enter the vulnerable
        // library path on the original multi-extension source.
        std::vector<RawHdu> raw;
        bool suspicious = false;
        bool incomplete = false;
        try {
            raw = scanRawHdus(source, incomplete);
            suspicious = incomplete;
            for (const auto& hdu : raw)
                suspicious |= hdu.bad_compression_header || hdu.repair_header;
        } catch (const std::exception& e) {
            // Do not hand an unparseable source to CFITSIO's full-HDU scan:
            // it may terminate the MPI rank with an uncaught native signal.
            throw std::runtime_error(std::string("raw FITS structural scan failed: ") + e.what());
        }
        std::vector<PlannedImage> plans;
        int invalid = incomplete ? 1 : 0;
        if (!suspicious) {
            int status = 0;
            fits_open_file(&input, source.string().c_str(), READONLY, &status);
            if (status != 0) throw std::runtime_error("cannot open source: " + fitsDiagnostic(status));
            try {
                plans = planImages(input, source, kind, final_directory, staging_directory, invalid);
            } catch (const std::exception& e) {
                if (raw.empty()) throw;
                std::cerr << "[init] " << source << " entering HDU salvage mode: " << e.what() << '\n';
                closeFits(input);
                invalid = incomplete ? 1 : 0; // Avoid double-counting fallback rescans.
                plans = rawPlans(source, kind, final_directory, staging_directory, raw, invalid);
            }
        } else {
            std::cerr << "[init] " << source << " malformed or dangerous FITS header detected; "
                      << "salvaging other HDUs\n";
            plans = rawPlans(source, kind, final_directory, staging_directory, raw, invalid);
        }
        result.skipped_hdus = invalid;
        if (plans.empty()) throw std::runtime_error("no safe extractable image HDUs");

        bool any_existing = false;
        for (const auto& plan : plans) any_existing |= fs::exists(plan.final_path);
        if (any_existing && policy == ExistingPolicy::Fail)
            throw std::runtime_error("output already exists: " + final_directory.string());

        std::vector<PlannedImage> pending;
        for (const auto& plan : plans) {
            if (fs::exists(plan.final_path) && policy == ExistingPolicy::Resume) {
                std::string error;
                if (validateOutput(plan, error)) {
                    result.output_paths.push_back(plan.final_path);
                    continue;
                }
                std::cerr << "[init] invalid resumed file " << plan.final_path
                          << ": " << error << "; re-extracting\n";
            }
            pending.push_back(plan);
        }
        if (pending.empty()) {
            closeFits(input);
            result.success = !result.output_paths.empty();
            result.resumed = result.success;
            return result;
        }
        fs::create_directories(staging_directory);
        fs::create_directories(final_directory);
        std::vector<PlannedImage> successful;
        for (const auto& plan : pending) {
            try {
                extractOne(source, plan, kind, staging_directory, input);
                successful.push_back(plan);
            } catch (const std::exception& e) {
                ++result.skipped_hdus;
                std::cerr << "[init] " << source << " HDU=" << plan.hdu_number
                          << " SKIP: " << e.what() << '\n';
                fits_clear_errmsg();
                std::error_code ec;
                fs::remove(plan.staged_path, ec);
            }
        }
        closeFits(input);
        const auto committed = commitImages(successful, source, result.skipped_hdus);
        for (const auto& plan : committed) result.output_paths.push_back(plan.final_path);
        result.success = !result.output_paths.empty();
        if (!result.success) result.error = "all planned HDUs failed extraction";
        std::error_code ec;
        fs::remove_all(staging_directory, ec);
        return result;
    } catch (const std::exception& e) {
        closeFits(input);
        result.error = e.what();
        std::error_code ec;
        fs::remove_all(staging_directory, ec);
        return result;
    }
}

} // namespace fqinit
