#include "process_main/PreProcess.hpp"
#include "process_main/PreProcessLegacy.hpp"
#include "process_main/ProcessMainState.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "pathconfig.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/Astrometry.hpp"
#include "process_main/ImageProcessing.hpp"
#include "general/NumericalRecipes.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace PreProcess {

    namespace {
        // ==========================================
        // Function: Estimate stripe median and width from valid blocks
        // Method: Preserve the F77 1000-sample robust estimator while drawing only from
        //         blocks that contain at least one usable pixel.
        // ==========================================
        void getStripeSigMed(const std::vector<float>& map,
                             const std::vector<unsigned char>& valid,
                             float& sig, float& med) {
            std::vector<float> values;
            values.reserve(map.size());
            for (size_t i = 0; i < map.size(); ++i) {
                if (valid[i] != 0U) values.push_back(map[i]);
            }

            if (values.empty()) {
                sig = -1.0f;
                med = 0.0f;
                return;
            }

            constexpr int sample_count = 1000;
            std::vector<float> samples(sample_count);
            for (int i = 0; i < sample_count; ++i) {
                size_t index = static_cast<size_t>(
                    NumericalRecipes::ran1() * static_cast<double>(values.size()));
                if (index >= values.size()) index = values.size() - 1;
                samples[static_cast<size_t>(i)] = values[index];
            }

            std::sort(samples.begin(), samples.end());
            sig = 0.5f * (samples[5 * sample_count / 6 - 1]
                          - samples[sample_count / 6 - 1]);
            med = samples[sample_count / 2 - 1];
        }
    }

    // Stage 1 driver
    void preProcess(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);
        
        for (const auto& image_file : image_files) {
            int cid = UniversalUtils::getChipId(image_file);

            chipPreProcess(image_file, dir_output, cid);
        }
    }

    // ==========================================
    // Function: Individual chip preprocessing
    // Method: Match the Fortran Stage 1 flow while keeping diagnostics side-effect free.
    // ==========================================
    void chipPreProcess(const std::string& imageFile, const std::string& dirOutput, int cid) {
        int proc_error = 0;
        int nx = 0, ny = 0;
        std::vector<float> array;
        WCSParams wcs;
        std::string prefix = UniversalUtils::getPrefix(imageFile);

        const bool image_read_ok = FitsIO::readImagePara(imageFile, nx, ny, array, wcs);
        const bool source_image_unreadable =
            !image_read_ok || (!array.empty() && array[0] < -99990.0f);
        if (source_image_unreadable) {
            std::cerr << "Error reading image parameters of: " << imageFile << std::endl;
            proc_error = 1;
            nx = 1;
            ny = 1;
            array.assign(1, -99999.0f);
        }

        // ==========================================
        // Function: Load the Stage-1 DQ mask before background and sigma fitting
        // Method: Keep DQ pixels in the processing-time weight mask while preserving their
        //         science values until the final norm.fits invalid-pixel serialization.
        // ==========================================
        std::vector<float> dqmask;
        if (proc_error == 0) {
            const std::string prefix_e = UniversalUtils::getPrefixExpo(imageFile);
            const std::string local_mask_file = dirOutput + "/dqmask/" + prefix_e + "/"
                                              + prefix_e + "_" + std::to_string(cid) + ".fits";
            int dnx = 0;
            int dny = 0;
            if (!FitsIO::readImage(local_mask_file, dnx, dny, dqmask)) {
                std::cerr << "Error / cant find mask file: " << local_mask_file << std::endl;
                proc_error = 1;
            } else if (dnx != nx || dny != ny) {
                std::cerr << "Error / wrong size of DQ file!" << std::endl;
                proc_error = 1;
            }
        }

        std::vector<int> weight(nx * ny, 1);
        std::vector<float> normap(nx * ny);

        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int idx = y * nx + x;
                if (array[idx] > LensingConfig::saturation_thresh) {
                    weight[idx] = 0;
                }
                if (proc_error == 0
                    && std::abs(dqmask[static_cast<size_t>(idx)]) > 1e-7f) {
                    weight[idx] = 0;
                }
                normap[idx] = array[idx];
            }
        }

        std::vector<double> bg_coeffs;
        std::vector<double> amplifier_bg_coeffs;
        std::vector<double> sig_coeffs;
        int nxc = nx / 2;
        double aa = 0.0, bb = 0.0, cc = 0.0;

        const auto collectBackground = [&](int x_start, int x_end) {
            amplifier_bg_coeffs.clear();
            setBackground(x_start, x_end, 0, ny, nx, ny, normap,
                          LensingConfig::blocksize, LensingConfig::nct, LensingConfig::ncx,
                          amplifier_bg_coeffs, proc_error);
            if (proc_error == 0) {
                bg_coeffs.insert(bg_coeffs.end(), amplifier_bg_coeffs.begin(),
                                 amplifier_bg_coeffs.end());
            }
        };

        if (LensingConfig::CCD_split == 2) {
            collectBackground(0, nxc);
            collectBackground(nxc, nx);
            
            setSig(0, nxc, 0, ny, nx, ny, normap, aa, bb, cc, proc_error);
            if (proc_error == 0) {
                sig_coeffs.push_back(aa);
                sig_coeffs.push_back(bb);
                sig_coeffs.push_back(cc);
            }
            if (proc_error == 0) {
                setSig(nxc, nx, 0, ny, nx, ny, normap, aa, bb, cc, proc_error);
            }
            if (proc_error == 0) {
                sig_coeffs.push_back(aa);
                sig_coeffs.push_back(bb);
                sig_coeffs.push_back(cc);
            }
        } else {
            collectBackground(0, nx);
            setSig(0, nx, 0, ny, nx, ny, normap, aa, bb, cc, proc_error);
            if (proc_error == 0) {
                sig_coeffs.push_back(aa);
                sig_coeffs.push_back(bb);
                sig_coeffs.push_back(cc);
            }
        }

        std::string astroFilename = OutputLayout::chipPath(
            dirOutput, "astrometry/dat_Astro", prefix, "_astro.dat");

        // ==========================================
        // Logic: Load the fixed repartitioned Gaia astrometry catalog layout
        // Method: Preserve the existing error sentinel, otherwise accumulate the
        //         configured one-degree candidates rooted at ASTROMETRY_CAT.
        // ==========================================
        if (proc_error != 0) {
            Astrometry::genAstrometryData(
                "", nx, ny, normap, weight, wcs, astroFilename, proc_error);
        } else {
            std::vector<std::string> catfiles = UniversalUtils::generateGalCatFileNames(
                LensingConfig::ASTROMETRY_CAT,
                wcs.crval,
                AstroCatConfig::ASTROMETRY_TILE_PREFIX);
            Astrometry::genAstrometryDataMulti(
                catfiles, nx, ny, normap, weight, wcs, astroFilename, proc_error);
        }

        locateDefects(nx, ny, array, normap, weight, LensingConfig::area_max, LensingConfig::area_thresh, proc_error);
        mergeDefects(nx, ny, weight, normap, LensingConfig::area_max, LensingConfig::source_thresh, LensingConfig::area_thresh, proc_error);

        if (proc_error == 0) {
            // Serialize the final combined saturation, DQ, and detected-defect mask.
            for (int i = 0; i < nx * ny; ++i) {
                if (weight[static_cast<size_t>(i)] == 0) {
                    normap[static_cast<size_t>(i)] = -1000.0f;
                }
            }
        } else {
            std::fill(normap.begin(), normap.end(), -1000.0f);
        }
        if (proc_error == 0) {
            normap[0] = -1.0f;
        } else {
            normap[0] = 1.0f;
        }

        std::string normFilename = OutputLayout::chipPath(
            dirOutput, "stamps/Norm", prefix, "_norm.fits");
        const bool norm_write_ok = source_image_unreadable
            ? FitsIO::writeImage(normFilename, 1, 1, std::vector<float>{1.0f})
            : FitsIO::writeNormHDU(imageFile, normFilename, nx, ny, normap,
                                   bg_coeffs, sig_coeffs, LensingConfig::CCD_split,
                                   LensingConfig::nct);
        if (!norm_write_ok) {
            std::cerr << "Error writing normalized image: " << normFilename << std::endl;
            proc_error = 1;
        }

        // if (proc_error == 0) {
        //     std::cout << "Status of processing " << imageFile << ": OK." << std::endl;
        // } else {
        //     std::cout << "Status of processing " << imageFile << ": ERROR!" << std::endl;
        // }
    }

    // ==========================================
    // Function: Estimate and subtract one amplifier background
    // Method: Delegate to the frozen historical F77 estimator and publish its metadata.
    // ==========================================
    void setBackground(int x_start, int x_end, int y_start, int y_end, int nx, int ny,
                       std::vector<float>& image, int blocksize, int nct, int ncx,
                       std::vector<double>& bg_coeffs, int& ierror) {
        bg_coeffs.clear();
        if (ierror != 0) return;
        if (!PreProcessLegacy::setBackground(
                x_start, x_end, y_start, y_end, nx, ny, image,
                blocksize, nct, ncx, bg_coeffs)) {
            ierror = 1;
        }
    }

    // ==========================================
    // Function: Estimate and apply one amplifier noise-sigma plane
    // Method: Delegate to the frozen historical F77 random-triple estimator.
    // ==========================================
    void setSig(int x_start, int x_end, int y_start, int y_end, int nx, int ny,
                std::vector<float>& image, double& aa, double& bb, double& cc,
                int& ierror) {
        aa = 0.0;
        bb = 0.0;
        cc = 0.0;
        if (ierror != 0) return;
        if (!PreProcessLegacy::setSig(
                x_start, x_end, y_start, y_end, nx, ny, image,
                aa, bb, cc)) {
            ierror = 1;
        }
    }

    // ==========================================
    // Function: Locate detector defects in one preprocessed CCD
    // Method: Stop before margin indexing whenever an earlier Stage-1 error is active.
    // ==========================================
    void locateDefects(int nx, int ny, const std::vector<float>& array, std::vector<float>& normap,
                       std::vector<int>& weight, int area_max, int area_thresh, int& ierror) {
        if (ierror != 0) return;

        constexpr int margin = 10;
        constexpr double defect_halo_thresh = 1.0;
        constexpr int y_smooth = 200;
        constexpr int x_smooth = 100;

        for (int y = 0; y < ny; ++y) {
            int mid_start = nx / 2 - margin - 1;
            int mid_end = nx / 2 + margin;
            for (int x = mid_start; x < mid_end; ++x) {
                if (x >= 0 && x < nx) {
                    weight[y * nx + x] = 0;
                }
            }
            for (int x = 0; x < margin; ++x) {
                weight[y * nx + x] = 0;
            }
            for (int x = nx - margin; x < nx; ++x) {
                weight[y * nx + x] = 0;
            }
        }
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < margin; ++y) {
                weight[y * nx + x] = 0;
            }
            for (int y = ny - margin; y < ny; ++y) {
                weight[y * nx + x] = 0;
            }
        }

       std::vector<float> map(nx * ny);

       for (int i = 0; i < nx * ny; ++i) {
           map[i] = UniversalUtils::loga(array[i], 1);
       }

        ImageProcessing::removeContinuous(nx, ny, nx, ny, map, UniversalUtils::iden, 4);

        std::vector<float> diffx(nx * ny);
        std::vector<float> diffy(nx * ny);
        for (int y = 0; y < ny; ++y) {
            int next_y = (y + 1) % ny;
            for (int x = 0; x < nx; ++x) {
                int next_x = (x + 1) % nx;
                diffx[y * nx + x] = map[y * nx + x] - map[y * nx + next_x];
                diffy[y * nx + x] = map[y * nx + x] - map[next_y * nx + x];
            }
        }

       float sigx = 0.0f, medx = 0.0f;
       ImageProcessing::getSigMed(nx, ny, diffx, sigx, medx);
       for (int i = 0; i < nx * ny; ++i) {
           if (std::abs(diffx[i]) > 8.0f * sigx) {
               weight[i] = 0;
           }
       }

       float sigy = 0.0f, medy = 0.0f;
       ImageProcessing::getSigMed(nx, ny, diffy, sigy, medy);
       for (int i = 0; i < nx * ny; ++i) {
           if (std::abs(diffy[i]) > 8.0f * sigy) {
               weight[i] = 0;
           }
       }

       maskSourceRegions(nx, ny, weight, normap, area_max, defect_halo_thresh * 2.0, area_thresh);

       detectStripes(nx, ny, normap, weight, x_smooth, y_smooth);

       detectArtificialStripes(nx, ny, weight, diffx, diffy, sigx, sigy, medx, medy);

       for (int i = 0; i < nx * ny; ++i) {
           if (weight[i] > 1) {
               weight[i] = 1;
           }
       }

       detectStellarHalo(nx, ny, normap, weight, area_max, defect_halo_thresh);

       detectDent(nx, ny, normap, weight, area_max, defect_halo_thresh);
    }

    // ==========================================
    // Function: Merge defect-connected source regions
    // Method: Traverse connected components in F77 x-major/u-major order and mask oversized regions.
    // ==========================================
   void mergeDefects(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                     int area_max, double source_thresh, int area_thresh, int& ierror) {
       if (ierror == 1) return;

       std::vector<int> mark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] >= source_thresh && weight[i] == 1) {
                mark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (mark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    mark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (mark[nidx] == 1) {
                                    component.push_back(nidx);
                                    mark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(area_max)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (mark[nidx] > 1 || (mark[nidx] == 0 && weight[nidx] == 0)) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            mark[p] = area_max;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            mark[p] = nb;
                        }
                    }
               }
           }
       }
   }

   void detectArtificialStripes(int nx, int ny, std::vector<int>& weight,
                                 const std::vector<float>& diffx, const std::vector<float>& diffy,
                                 float sigx, float sigy, float medx, float medy) {
        std::vector<float> entropy(nx * ny);
        ImageProcessing::getEntropy(nx, ny, diffx, sigx, medx, 2, entropy);
        
        float sig = 0.0f, med = 0.0f;
        ImageProcessing::getSigMed(nx, ny, entropy, sig, med);
        for (int i = 0; i < nx * ny; ++i) {
            if (weight[i] == 1 && std::abs(entropy[i] - med) > 10.0f * sig) {
                weight[i] = 0;
            }
        }

        ImageProcessing::getEntropy(nx, ny, diffy, sigy, medy, 2, entropy);
        ImageProcessing::getSigMed(nx, ny, entropy, sig, med);
        for (int i = 0; i < nx * ny; ++i) {
            if (weight[i] == 1 && std::abs(entropy[i] - med) > 10.0f * sig) {
                weight[i] = 0;
            }
        }
    }

    // ==========================================
    // Function: Mark bright source regions before stripe and halo detection
    // Method: Match F77 mask_source_regions connected-component traversal and area_max handling.
    // ==========================================
    void maskSourceRegions(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                           int area_max, double source_thresh, int area_thresh) {
        std::vector<int> mark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] >= source_thresh && weight[i] == 1) {
                mark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (mark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    mark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (mark[nidx] == 1) {
                                    component.push_back(nidx);
                                    mark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(area_max)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (mark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            mark[p] = area_max;
                            weight[p] = 2;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        if (nb >= area_thresh) {
                            for (int p : component) {
                                mark[p] = nb;
                                weight[p] = 2;
                            }
                        } else {
                            for (int p : component) {
                                mark[p] = nb;
                            }
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect stripe outliers including partial edge blocks
    // Method: Use ceiling block counts and sum/sqrt(nvalid), so full, partial, and DQ-masked
    //         blocks share one noise scale.
    // ==========================================
    void detectStripes(int nx, int ny, const std::vector<float>& normap,
                       std::vector<int>& weight, int x_smooth, int y_smooth) {
        int numy = (ny + y_smooth - 1) / y_smooth;
        std::vector<float> ymap(nx * numy, 0.0f);
        std::vector<unsigned char> yvalid(nx * numy, 0U);
        for (int x = 0; x < nx; ++x) {
            for (int y_block = 0; y_block < numy; ++y_block) {
                int starty = y_block * y_smooth;
                int endy = std::min(starty + y_smooth, ny);
                float sum = 0.0f;
                int nvalid = 0;
                for (int y = starty; y < endy; ++y) {
                    if (weight[y * nx + x] == 1) {
                        sum += normap[y * nx + x];
                        ++nvalid;
                    }
                }
                if (nvalid > 0) {
                    size_t index = static_cast<size_t>(y_block * nx + x);
                    ymap[index] = sum / std::sqrt(static_cast<float>(nvalid));
                    yvalid[index] = 1U;
                }
            }
        }

        float sig_y = 0.0f, med_y = 0.0f;
        getStripeSigMed(ymap, yvalid, sig_y, med_y);
        for (int x = 0; x < nx; ++x) {
            for (int y_block = 0; y_block < numy; ++y_block) {
                size_t index = static_cast<size_t>(y_block * nx + x);
                if (yvalid[index] != 0U && sig_y > 0.0f
                    && std::abs(ymap[index] - med_y) > sig_y * 4.0f) {
                    int starty = y_block * y_smooth;
                    int endy = std::min(starty + y_smooth, ny);
                    for (int y = starty; y < endy; ++y) {
                        weight[y * nx + x] = 0;
                    }
                }
            }
        }

        int numx = (nx + x_smooth - 1) / x_smooth;
        std::vector<float> xmap(numx * ny, 0.0f);
        std::vector<unsigned char> xvalid(numx * ny, 0U);
        for (int y = 0; y < ny; ++y) {
            for (int x_block = 0; x_block < numx; ++x_block) {
                int startx = x_block * x_smooth;
                int endx = std::min(startx + x_smooth, nx);
                float sum = 0.0f;
                int nvalid = 0;
                for (int x = startx; x < endx; ++x) {
                    if (weight[y * nx + x] == 1) {
                        sum += normap[y * nx + x];
                        ++nvalid;
                    }
                }
                if (nvalid > 0) {
                    size_t index = static_cast<size_t>(y * numx + x_block);
                    xmap[index] = sum / std::sqrt(static_cast<float>(nvalid));
                    xvalid[index] = 1U;
                }
            }
        }

        float sig_x = 0.0f, med_x = 0.0f;
        getStripeSigMed(xmap, xvalid, sig_x, med_x);
        for (int y = 0; y < ny; ++y) {
            for (int x_block = 0; x_block < numx; ++x_block) {
                size_t index = static_cast<size_t>(y * numx + x_block);
                if (xvalid[index] != 0U && sig_x > 0.0f
                    && std::abs(xmap[index] - med_x) > sig_x * 4.0f) {
                    int startx = x_block * x_smooth;
                    int endx = std::min(startx + x_smooth, nx);
                    for (int x = startx; x < endx; ++x) {
                        weight[y * nx + x] = 0;
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect stellar halo affected regions
    // Method: Match F77 detect_stellar_halo by smoothing the full normalized map before threshold-connected component masking.
    // ==========================================
    void detectStellarHalo(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                           int npmax, double defect_halo_thresh) {
        std::vector<float> smoothed = normap;
        ImageProcessing::smoothImage55(nx, ny, smoothed, 1);

        std::vector<int> dmark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (smoothed[i] >= defect_halo_thresh) {
                dmark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (dmark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    dmark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (dmark[nidx] == 1) {
                                    component.push_back(nidx);
                                    dmark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(npmax)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (dmark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            dmark[p] = npmax;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            dmark[p] = nb;
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Detect large negative dent regions
    // Method: Match F77 detect_dent connected-component traversal and area_max masking.
    // ==========================================
    void detectDent(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                    int npmax, double defect_halo_thresh) {
        std::vector<int> dmark(nx * ny, 0);
        for (int i = 0; i < nx * ny; ++i) {
            if (normap[i] <= -defect_halo_thresh) {
                dmark[i] = 1;
            }
        }

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                int idx = y * nx + x;
                if (dmark[idx] == 1) {
                    std::vector<int> component;
                    component.push_back(idx);
                    dmark[idx] = -1;

                    bool toobig = false;
                    size_t head = 0;
                    while (head < component.size()) {
                        int curr = component[head++];
                        int cx = curr % nx;
                        int cy = curr / nx;

                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx_val = cx + dx;
                            if (nx_val < 0 || nx_val >= nx) continue;
                            for (int dy = -1; dy <= 1; ++dy) {
                                int ny_val = cy + dy;
                                if (ny_val < 0 || ny_val >= ny) continue;

                                int nidx = ny_val * nx + nx_val;
                                if (dmark[nidx] == 1) {
                                    component.push_back(nidx);
                                    dmark[nidx] = -1;
                                    if (component.size() == static_cast<size_t>(npmax)) {
                                        toobig = true;
                                        break;
                                    }
                                } else if (dmark[nidx] > 1) {
                                    toobig = true;
                                    break;
                                }
                            }
                            if (toobig) break;
                        }
                        if (toobig) break;
                    }

                    if (toobig) {
                        for (int p : component) {
                            dmark[p] = npmax;
                            weight[p] = 0;
                        }
                    } else {
                        int nb = static_cast<int>(component.size());
                        for (int p : component) {
                            dmark[p] = nb;
                        }
                    }
                }
            }
        }
    }

}
