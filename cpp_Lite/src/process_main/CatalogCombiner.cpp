#include "process_main/CatalogCombiner.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/MPIFailure.hpp"
#include "process_main/OutputFile.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/ExposureInfo.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <system_error>

namespace CatalogCombiner {

// ==========================================
// Function: Count every physical line in one Stage-9 input catalog
// Method: Use a fresh getline stream and distinguish clean EOF from an I/O failure.
// ==========================================
static std::size_t countCatalogLines(const std::string& filename,
                                     const std::string& role) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        MPIFailure::abortWorld(
            "open catalog for row-count preflight", role + "=" + filename);
    }

    std::size_t line_count = 0;
    std::string line;
    while (std::getline(input, line)) {
        ++line_count;
    }
    if (input.bad()) {
        MPIFailure::abortWorld(
            "read catalog for row-count preflight", role + "=" + filename);
    }
    return line_count;
}

// ==========================================
// Function: Determine the fixed number of paired Stage-9 data rows
// Method: Count shear then orig, retry both with fresh streams after one mismatch,
//         and preserve a one-line shear catalog as the header-only sentinel.
// ==========================================
static std::size_t determinePairedDataRows(
    const std::string& filename_shear,
    const std::string& filename_orig,
    const std::string& prefix) {
    const std::size_t shear_1 = countCatalogLines(filename_shear, "shear");
    if (shear_1 == 0) {
        MPIFailure::abortWorld(
            "preflight Stage 7 shear catalog",
            "shear catalog contains no header prefix=" + prefix
                + " shear=" + filename_shear);
    }
    if (shear_1 == 1) {
        return 0;
    }

    const std::size_t orig_1 = countCatalogLines(filename_orig, "orig");
    if (shear_1 == orig_1) {
        return shear_1 - 1;
    }

    const std::size_t shear_2 = countCatalogLines(filename_shear, "shear");
    if (shear_2 == 0) {
        MPIFailure::abortWorld(
            "preflight Stage 7 shear catalog",
            "shear catalog contains no header prefix=" + prefix
                + " shear=" + filename_shear);
    }
    if (shear_2 == 1) {
        return 0;
    }

    const std::size_t orig_2 = countCatalogLines(filename_orig, "orig");
    if (shear_2 != orig_2) {
        std::ostringstream detail;
        detail << "prefix=" << prefix
               << " attempt1_shear_lines=" << shear_1
               << " attempt1_orig_lines=" << orig_1
               << " attempt2_shear_lines=" << shear_2
               << " attempt2_orig_lines=" << orig_2
               << " shear=" << filename_shear
               << " orig=" << filename_orig;
        MPIFailure::abortWorld(
            "combine catalog row-count preflight", detail.str());
    }

    std::cout << "CATALOG_ROWCOUNT_RECOVERED"
              << " prefix=" << prefix
              << " attempt1_shear_lines=" << shear_1
              << " attempt1_orig_lines=" << orig_1
              << " attempt2_shear_lines=" << shear_2
              << " attempt2_orig_lines=" << orig_2 << std::endl;
    return shear_2 - 1;
}

// ==========================================
// Function: Remove trailing ASCII whitespace from one catalog line
// Method: Pop only the established space, CR, LF, and tab suffix characters.
// ==========================================
static inline std::string trimRight(std::string str) {
    while (!str.empty() && (str.back() == ' ' || str.back() == '\r' || str.back() == '\n' || str.back() == '\t')) {
        str.pop_back();
    }
    return str;
}

// ==========================================
// Function: Combine one exposure's chip catalogs into the final result catalog
// Method: Preflight paired physical row counts with one fresh retry, preserve
//         header-only shear sentinels, and consume the matched rows in a fixed loop.
// ==========================================
void combineExpoCatalog(int nchip,
                        const std::vector<std::string>& imageFiles,
                        const std::string& dirOutput,
                        int expo_index,
                        float chi2) {
    if (nchip <= 0 || imageFiles.empty()) {
        MPIFailure::abortWorld(
            "combine exposure catalog", "exposure contains no chip paths");
    }

    std::string prefix_expo = UniversalUtils::getPrefixExpo(imageFiles[0]);
    std::string out_filename = dirOutput + "/result/" + prefix_expo + "_all.cat";

    std::error_code remove_error;
    std::filesystem::remove(out_filename, remove_error);
    if (remove_error) {
        MPIFailure::abortWorld(
            "remove stale combined catalog",
            out_filename + ": " + remove_error.message());
    }

    MainIO::OutputFile fout20;
    bool output_opened = false;

    int n = 0;
    int m = 0;
    int num_cols = LensingConfig::iparity + 1; // 24 columns

    std::string last_prefix;

    for (int ichip = 0; ichip < nchip; ++ichip) {
        const Universalblock::NormStatus norm_status =
            Universalblock::checkNorm(imageFiles[ichip], dirOutput);
        if (norm_status == Universalblock::NormStatus::Invalid) {
            continue;
        }
        if (norm_status != Universalblock::NormStatus::Valid) {
            MPIFailure::abortWorld(
                "check Stage 1 norm before catalog combination",
                Universalblock::normErrorDetail(
                    norm_status, imageFiles[ichip], dirOutput));
        }

        int chip_index = UniversalUtils::getChipId(imageFiles[ichip]);
        std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
        last_prefix = prefix;

        std::string filename_shear = OutputLayout::chipPath(
            dirOutput, "stamps/dat_Shear", prefix, "_shear.dat");
        std::string filename_orig = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
        const std::size_t paired_data_rows = determinePairedDataRows(
            filename_shear, filename_orig, prefix);
        if (paired_data_rows == 0) {
            continue;
        }

        std::ifstream fin10(filename_shear);
        if (!fin10.is_open()) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }

        std::string cat_list1;
        if (!std::getline(fin10, cat_list1)) {
            MPIFailure::abortWorld(
                "read Stage 7 shear catalog header", filename_shear);
        }
        cat_list1 = trimRight(cat_list1);
        if (cat_list1.empty()) {
            MPIFailure::abortWorld(
                "parse Stage 7 shear catalog header", filename_shear);
        }

        std::ifstream fin15(filename_orig);
        if (!fin15.is_open()) {
            MPIFailure::abortWorld("read external source catalog", filename_orig);
        }

        std::string cat_list2;
        if (!std::getline(fin15, cat_list2)) {
            MPIFailure::abortWorld(
                "read external source catalog header", filename_orig);
        }
        cat_list2 = trimRight(cat_list2);
        if (cat_list2.empty()) {
            MPIFailure::abortWorld(
                "parse external source catalog header", filename_orig);
        }

        if (chi2 > LensingConfig::chi2_thresh) {
            std::cout << prefix << " contains no valid sources!" << std::endl;
            return;
        }

        if (!output_opened) {
            fout20.open(out_filename);
            fout20 << std::setprecision(10);
            fout20 << cat_list2 << " EXPO_NUM ccD_NUM "
                   << cat_list1 << " Chi2\n";
            output_opened = true;
        }

        for (std::size_t row_index = 0;
             row_index < paired_data_rows; ++row_index) {
            std::string line10;
            std::string cat_content;
            const bool shear_ok = static_cast<bool>(std::getline(fin10, line10));
            const bool orig_ok = static_cast<bool>(std::getline(fin15, cat_content));
            if (!shear_ok || !orig_ok) {
                std::ostringstream detail;
                detail << "prefix=" << prefix
                       << " row_index_zero_based=" << row_index
                       << " row_index_one_based=" << (row_index + 1)
                       << " expected_data_rows=" << paired_data_rows
                       << " shear_read_ok=" << shear_ok
                       << " orig_read_ok=" << orig_ok
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "read fixed paired catalog row", detail.str());
            }

            cat_content = trimRight(cat_content);
            if (cat_content.empty()) {
                std::ostringstream detail;
                detail << "empty external catalog data row prefix=" << prefix
                       << " pair_index=" << (row_index + 1)
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired external catalog row", detail.str());
            }

            if (line10.empty()) {
                std::ostringstream detail;
                detail << "empty shear data row prefix=" << prefix
                       << " pair_index=" << (row_index + 1)
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired Stage 7 shear row", detail.str());
            }
            std::stringstream ss(line10);
            std::vector<float> cat(num_cols);
            bool read_ok = true;
            for (int u = 0; u < num_cols; ++u) {
                if (!(ss >> cat[u])) {
                    read_ok = false;
                    break;
                }
            }
            if (!read_ok) {
                std::ostringstream detail;
                detail << "incomplete shear data row prefix=" << prefix
                       << " pair_index=" << (row_index + 1)
                       << " expected_columns=" << num_cols
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired Stage 7 shear row", detail.str());
            }

            // A parsed shear row and its original catalog row form one pair.
            // Consume both before applying scientific rejection so the streams
            // remain aligned even when the pair is omitted from final output.
            if (cat[LensingConfig::i_imax] >= LensingConfig::ns || cat[LensingConfig::i_jmax] >= LensingConfig::ns) {
                m++;
                continue;
            }
            if (cat[0] < -900.0f) {
                m++;
                continue;
            }

            n++;
            double g1c = 0.0;
            double g2c = 0.0;
            // double g1c = static_cast<double>(cat[LensingConfig::igf1]) + LensingConfig::g1_c;
            // double g2c = static_cast<double>(cat[LensingConfig::igf2]) + LensingConfig::g2_c;
            cat[LensingConfig::ig1] = static_cast<float>(cat[LensingConfig::ig1] - g1c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih1] + g2c * cat[LensingConfig::ih2]);
            cat[LensingConfig::ig2] = static_cast<float>(cat[LensingConfig::ig2] - g2c * cat[LensingConfig::ide] + g1c * cat[LensingConfig::ih2] - g2c * cat[LensingConfig::ih1]);

            fout20 << cat_content << " " << expo_index << " " << chip_index;
            for (int u = 0; u < num_cols; ++u) {
                fout20 << " " << cat[u];
            }
            fout20 << " " << chi2 << "\n";
        }
    }

    std::cout << (last_prefix.empty() ? prefix_expo : last_prefix)
              << " " << n << " " << m << std::endl;
    if (output_opened) {
        fout20.close();
    }
}

// ==========================================
// Function: Run Stage 9 for one exposure
// Method: Resolve chip paths, obtain reduced Stage-8 chi2, and invoke the combiner.
// ==========================================
void procComb(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    float chi2 = 0.0f;
    if (ExposureInfo::state.parameters.size() >= static_cast<size_t>(iexpo) * 6) {
        chi2 = ExposureInfo::state.parameters[(iexpo - 1) * 6 + 2]; // 3rd element in Fortran, index 2
    }

    combineExpoCatalog(
        static_cast<int>(image_files.size()), image_files, dir_output, iexpo, chi2);
}

} // namespace CatalogCombiner
