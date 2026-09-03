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

static inline std::string trimRight(std::string str) {
    while (!str.empty() && (str.back() == ' ' || str.back() == '\r' || str.back() == '\n' || str.back() == '\t')) {
        str.pop_back();
    }
    return str;
}

// ==========================================
// Function: Combine one exposure's chip catalogs into the final result catalog
// Method: Remove stale output, preserve Stage-7 header-only sentinels, consume
//         live shear/orig records in lockstep, and lazily open the final catalog.
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

        std::string line10;
        if (!std::getline(fin10, line10)) {
            continue;
        }

        std::string filename_orig = OutputLayout::chipPath(
            dirOutput, "stamps/cat_Orig", prefix, "_orig.cat");
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

        std::size_t pair_index = 0;
        do {
            ++pair_index;
            if (line10.empty()) {
                std::ostringstream detail;
                detail << "empty shear data row prefix=" << prefix
                       << " pair_index=" << pair_index
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
                       << " pair_index=" << pair_index
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired Stage 7 shear row", detail.str());
            }

            std::string cat_content;
            if (!std::getline(fin15, cat_content)) {
                std::ostringstream detail;
                detail << "external catalog ended before shear catalog"
                       << " prefix=" << prefix
                       << " pair_index=" << pair_index
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld("combine paired catalog", detail.str());
            }
            cat_content = trimRight(cat_content);
            if (cat_content.empty()) {
                std::ostringstream detail;
                detail << "empty external catalog data row prefix=" << prefix
                       << " pair_index=" << pair_index
                       << " shear=" << filename_shear
                       << " orig=" << filename_orig;
                MPIFailure::abortWorld(
                    "parse paired external catalog row", detail.str());
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
        } while (std::getline(fin10, line10));

        if (fin10.bad()) {
            MPIFailure::abortWorld("read Stage 7 shear catalog", filename_shear);
        }
        std::string extra_orig;
        if (std::getline(fin15, extra_orig)) {
            std::ostringstream detail;
            detail << "external catalog contains unpaired trailing row"
                   << " prefix=" << prefix
                   << " pair_index=" << (pair_index + 1)
                   << " shear=" << filename_shear
                   << " orig=" << filename_orig;
            MPIFailure::abortWorld("combine paired catalog", detail.str());
        }
        if (fin15.bad()) {
            MPIFailure::abortWorld(
                "read external source catalog", filename_orig);
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
