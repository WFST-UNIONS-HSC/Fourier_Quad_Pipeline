#include "process_fd/ShearCatalogReader.hpp"
#include "FDConfig.hpp"
#include "general/NumericalRecipes.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

namespace fc = FDConfig;

// ==========================================
// Function: readExposure
// Method: Group raw valid rows by external-catalog coordinates, uniformly
//         retain at most MAX_DUP measurements, then apply the existing FD cuts
//         and append accepted measurements to the shared FDData arrays.
// ==========================================
void ShearCatalogReader::readExposure(int iexpo, FDData& data,
                                      const std::vector<std::string>& expo_files,
                                      int rank) {
    if (iexpo < 1 || iexpo > static_cast<int>(expo_files.size())) {
        if (rank == 0)
            std::cerr << "Invalid exposure index: " << iexpo << std::endl;
        return;
    }

    // ==========================================
    // Logic: Derive the accepted catalog footprint from the configured CCD geometry
    // Method: Apply one symmetric edge width instead of fixed survey-specific maxima.
    // ==========================================
    const int chip_xmin = fc::chip_mask_edge;
    const int chip_xmax = LensingConfig::chipnx - fc::chip_mask_edge;
    const int chip_ymin = fc::chip_mask_edge;
    const int chip_ymax = LensingConfig::chipny - fc::chip_mask_edge;

    const std::string& filename = expo_files[iexpo - 1];
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << filename << " does not exist!!" << std::endl;
        return;
    }

    // Skip header line
    std::string header;
    std::getline(file, header);

    // Duplicate-detection buffer (stores up to MAX_DUP rows per source)
    std::vector<std::vector<float>> dup_buf(fc::MAX_DUP,
                                            std::vector<float>(fc::ICHI2));
    int ndup = 0;
    float last_dec = -999.0, last_ra = -999.0;
    bool have_group = false;
    constexpr float multisam_thrsh = 1e-7f;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<float> item(fc::ICHI2);
        bool parse_ok = true;
        for (int u = 0; u < fc::ICHI2; ++u) {
            if (!(iss >> item[u])) { parse_ok = false; break; }
        }
        if (!parse_ok) continue;

        // NaN / Inf check
        bool undefined = false;
        for (int u = 0; u < fc::ICHI2; ++u) {
            if (std::isnan(item[u])) { undefined = true; break; }
            if (std::fabs(item[u]) > 1.0e30) { undefined = true; break; }
        }
        if (undefined) continue;

        // Duplicate grouping intentionally precedes every scientific cut.
        const float current_dec = item[fc::col_cdec];
        const float current_ra = item[fc::col_cra];
        const bool new_galaxy =
            have_group
            && (std::fabs(current_dec - last_dec) > multisam_thrsh
                || std::fabs(current_ra - last_ra) > multisam_thrsh);

        // Flush the sampled measurements from the completed galaxy group.
        if (new_galaxy && ndup > 0) {
            const int selected_count =
                ndup > fc::MAX_DUP ? fc::MAX_DUP : ndup;
            for (int i = 0; i < selected_count; ++i) {
                const std::vector<float>& row = dup_buf[i];

                int ix = static_cast<int>(row[fc::col_pixx]);
                int iy = static_cast<int>(row[fc::col_pixy]);
                int ccd_val = static_cast<int>(std::lround(row[fc::col_ccd]));
                bool bad_ccd = false;
                for (int j = 0; j < fc::n_bad_ccds; ++j) {
                    if (ccd_val == fc::bad_ccds[j]) {
                        bad_ccd = true;
                        break;
                    }
                }
                if (bad_ccd) continue;
                if (ix < chip_xmin || ix > chip_xmax ||
                    iy < chip_ymin || iy > chip_ymax) continue;
                if (fc::ft_cut >= 0.0 &&
                    std::fabs(row[fc::col_flags_ft] - fc::ft_cut) > 1e-3) continue;
                if (fc::fg_cut >= 0.0 &&
                    std::fabs(row[fc::col_flags_fg] - fc::fg_cut) > 1e-3) continue;
                if (fc::gold_cut >= 0.0 &&
                    std::fabs(row[fc::col_flags_gold] - fc::gold_cut) > 1e-3) continue;
                if (fc::ext_cut >= 0.0 &&
                    std::fabs(row[fc::col_ext_mash] - fc::ext_cut) > 1e-3) continue;
                if (row[fc::col_SNR_F] < fc::snrfcut) continue;

                float snr = row[fc::col_h_flux] /
                            std::sqrt(row[fc::col_h_area]);
                if (fc::snrlow > 0.0 && snr < fc::snrlow) continue;
                if (fc::snrhigh > 0.0 && snr > fc::snrhigh) continue;
                if (fc::r_half_thresh > 0.0) {
                    float r_half =
                        std::sqrt(row[fc::col_h_area] / LensingConfig::pi) *
                        LensingConfig::pixel_size * 2.0;
                    if (r_half <= fc::r_half_thresh * row[fc::col_PSF]) continue;
                }
                if (row[fc::col_mag_i] < 10.0 || row[fc::col_mag_i] > 30.0) continue;
                if (row[fc::col_polychi2] > fc::psf_chi2_mltp) continue;
                if (row[fc::col_star] < fc::starcut) continue;
                if (!fc::FD_PER_EXPOSURE_STAR_BAR &&
                    row[fc::col_chi2] > fc::chi2_thresh) continue;
                if (row[fc::col_flag] <= fc::flagcut) continue;
                if (row[fc::col_imax] >= fc::imaxcut) continue;
                if (row[fc::col_jmax] >= fc::jmaxcut) continue;
                if (row[fc::col_zp] <= fc::zplow) continue;
                if (row[fc::col_zp] >= fc::zphigh) continue;
                if (std::fabs(row[fc::col_gf1]) > 0.0015) continue;
                if (std::fabs(row[fc::col_gf2]) > 0.0015) continue;

                // ==========================================
                // Logic: Preserve the source's original exposure identity
                // Method: Read the 1-based catalog value and reject IDs that
                //         cannot address an exposure-sized vector.
                // ==========================================
                int source_expo = 0;
                if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                    const double exposure_value = row[fc::col_expo];
                    if (!std::isfinite(exposure_value)
                        || exposure_value < 1.0
                        || exposure_value > static_cast<double>(std::numeric_limits<int>::max())) {
                        continue;
                    }
                    source_expo = static_cast<int>(std::lround(exposure_value));
                }

                int idx = data.ng;
                if (idx >= fc::nmax_per_core) {
                    std::cerr << "nmax_per_core is too small!" << std::endl;
                    file.close();
                    return;
                }
                data.x1[idx]  = dup_buf[i][fc::col_gf1];
                data.y1[idx]  = dup_buf[i][fc::col_g1];
                data.de1[idx] = dup_buf[i][fc::col_de] - dup_buf[i][fc::col_h1];
                data.x2[idx]  = dup_buf[i][fc::col_gf2];
                data.y2[idx]  = dup_buf[i][fc::col_g2];
                data.de2[idx] = dup_buf[i][fc::col_de] + dup_buf[i][fc::col_h1];

                // Jackknife weight: ww = ±1/sqrt((g1/de)² + (g2/de)²)
                {
                    float de_val = dup_buf[i][fc::col_de];
                    float y1j = dup_buf[i][fc::col_g1] / de_val;
                    float y2j = dup_buf[i][fc::col_g2] / de_val;
                    float gmag = std::sqrt(y1j * y1j + y2j * y2j);
                    float sign = (de_val >= 0.0) ? 1.0 : -1.0;
                    data.ww[idx] = (gmag > 0.0) ? sign / gmag : 0.0;
                }

                data.magr[idx]    = dup_buf[i][fc::col_mag_r];
                data.magg[idx]    = dup_buf[i][fc::col_mag_g];
                data.magi[idx]    = dup_buf[i][fc::col_mag_i];
                data.sizerel[idx] =
                    ((std::sqrt(dup_buf[i][fc::col_h_area] / LensingConfig::pi) * LensingConfig::pixel_size * 2.0)
                     - dup_buf[i][fc::col_PSF]) / dup_buf[i][fc::col_PSF];
                data.src_snr[idx] = dup_buf[i][fc::col_h_flux] /
                                     std::sqrt(dup_buf[i][fc::col_h_area]);
                data.rra[idx]  = dup_buf[i][fc::col_ra];
                data.ddec[idx] = dup_buf[i][fc::col_dec];

                if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                    data.iexpo[idx] = source_expo;
                    data.snrf[idx]  = dup_buf[i][fc::col_SNR_F];
                }
                data.ng++;
            }
            ndup = 0;
        }

        if (!have_group || new_galaxy) {
            last_dec = current_dec;
            last_ra = current_ra;
            have_group = true;
        }

        // ==========================================
        // Logic: Uniformly retain MAX_DUP measurements from the raw group
        // Method: Apply reservoir sampling before any scientific cut, using
        //         the per-rank ran1 stream initialized once by main.cpp.
        // ==========================================
        if (ndup < fc::MAX_DUP) {
            dup_buf[ndup] = item;
        } else {
            const int reservoir_index = static_cast<int>(
                NumericalRecipes::ran1() * static_cast<double>(ndup + 1));
            if (reservoir_index < fc::MAX_DUP) {
                dup_buf[reservoir_index] = item;
            }
        }
        ++ndup;
    }

    // Flush remaining duplicates
    if (ndup > 0) {
        const int selected_count =
            ndup > fc::MAX_DUP ? fc::MAX_DUP : ndup;
        for (int i = 0; i < selected_count; ++i) {
            const std::vector<float>& row = dup_buf[i];

            int ix = static_cast<int>(row[fc::col_pixx]);
            int iy = static_cast<int>(row[fc::col_pixy]);
            int ccd_val = static_cast<int>(std::lround(row[fc::col_ccd]));
            bool bad_ccd = false;
            for (int j = 0; j < fc::n_bad_ccds; ++j) {
                if (ccd_val == fc::bad_ccds[j]) {
                    bad_ccd = true;
                    break;
                }
            }
            if (bad_ccd) continue;
            if (ix < chip_xmin || ix > chip_xmax ||
                iy < chip_ymin || iy > chip_ymax) continue;
            if (fc::ft_cut >= 0.0 &&
                std::fabs(row[fc::col_flags_ft] - fc::ft_cut) > 1e-3) continue;
            if (fc::fg_cut >= 0.0 &&
                std::fabs(row[fc::col_flags_fg] - fc::fg_cut) > 1e-3) continue;
            if (fc::gold_cut >= 0.0 &&
                std::fabs(row[fc::col_flags_gold] - fc::gold_cut) > 1e-3) continue;
            if (fc::ext_cut >= 0.0 &&
                std::fabs(row[fc::col_ext_mash] - fc::ext_cut) > 1e-3) continue;
            if (row[fc::col_SNR_F] < fc::snrfcut) continue;

            float snr = row[fc::col_h_flux] /
                        std::sqrt(row[fc::col_h_area]);
            if (fc::snrlow > 0.0 && snr < fc::snrlow) continue;
            if (fc::snrhigh > 0.0 && snr > fc::snrhigh) continue;
            if (fc::r_half_thresh > 0.0) {
                float r_half =
                    std::sqrt(row[fc::col_h_area] / LensingConfig::pi) *
                    LensingConfig::pixel_size * 2.0;
                if (r_half <= fc::r_half_thresh * row[fc::col_PSF]) continue;
            }
            if (row[fc::col_mag_i] < 10.0 || row[fc::col_mag_i] > 30.0) continue;
            if (row[fc::col_polychi2] > fc::psf_chi2_mltp) continue;
            if (row[fc::col_star] < fc::starcut) continue;
            if (!fc::FD_PER_EXPOSURE_STAR_BAR &&
                row[fc::col_chi2] > fc::chi2_thresh) continue;
            if (row[fc::col_flag] <= fc::flagcut) continue;
            if (row[fc::col_imax] >= fc::imaxcut) continue;
            if (row[fc::col_jmax] >= fc::jmaxcut) continue;
            if (row[fc::col_zp] <= fc::zplow) continue;
            if (row[fc::col_zp] >= fc::zphigh) continue;
            if (std::fabs(row[fc::col_gf1]) > 0.0015) continue;
            if (std::fabs(row[fc::col_gf2]) > 0.0015) continue;

            // ==========================================
            // Logic: Preserve the source's original exposure identity
            // Method: Read the 1-based catalog value and reject IDs that
            //         cannot address an exposure-sized vector.
            // ==========================================
            int source_expo = 0;
            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                const double exposure_value = row[fc::col_expo];
                if (!std::isfinite(exposure_value)
                    || exposure_value < 1.0
                    || exposure_value > static_cast<double>(std::numeric_limits<int>::max())) {
                    continue;
                }
                source_expo = static_cast<int>(std::lround(exposure_value));
            }

            int idx = data.ng;
            if (idx >= fc::nmax_per_core) {
                std::cerr << "nmax_per_core is too small!" << std::endl;
                break;
            }
            data.x1[idx]  = dup_buf[i][fc::col_gf1];
            data.y1[idx]  = dup_buf[i][fc::col_g1];
            data.de1[idx] = dup_buf[i][fc::col_de] - dup_buf[i][fc::col_h1];
            data.x2[idx]  = dup_buf[i][fc::col_gf2];
            data.y2[idx]  = dup_buf[i][fc::col_g2];
            data.de2[idx] = dup_buf[i][fc::col_de] + dup_buf[i][fc::col_h1];

            // Jackknife weight: ww = ±1/sqrt((g1/de)² + (g2/de)²)
            {
                float de_val = dup_buf[i][fc::col_de];
                float y1j = dup_buf[i][fc::col_g1] / de_val;
                float y2j = dup_buf[i][fc::col_g2] / de_val;
                float gmag = std::sqrt(y1j * y1j + y2j * y2j);
                float sign = (de_val >= 0.0) ? 1.0 : -1.0;
                data.ww[idx] = (gmag > 0.0) ? sign / gmag : 0.0;
            }

            data.magr[idx]    = dup_buf[i][fc::col_mag_r];
            data.magg[idx]    = dup_buf[i][fc::col_mag_g];
            data.magi[idx]    = dup_buf[i][fc::col_mag_i];
            data.sizerel[idx] =
                ((std::sqrt(dup_buf[i][fc::col_h_area] / LensingConfig::pi) * LensingConfig::pixel_size * 2.0)
                 - dup_buf[i][fc::col_PSF]) / dup_buf[i][fc::col_PSF];
            data.src_snr[idx] = dup_buf[i][fc::col_h_flux] /
                                 std::sqrt(dup_buf[i][fc::col_h_area]);
            data.rra[idx]  = dup_buf[i][fc::col_ra];
            data.ddec[idx] = dup_buf[i][fc::col_dec];

            if (fc::FD_PER_EXPOSURE_STAR_BAR) {
                data.iexpo[idx] = source_expo;
                data.snrf[idx]  = dup_buf[i][fc::col_SNR_F];
            }
            data.ng++;
        }
    }
    file.close();
}
