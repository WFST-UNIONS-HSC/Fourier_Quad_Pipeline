#include "process_main/PSFModel.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/PSFModelState.hpp"
#include "process_main/PSFCandidateQuality.hpp"
#include "process_main/PSFStarSelection.hpp"
#include "process_main/OutputFile.hpp"
#include "process_main/MPIFailure.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/Astrometry.hpp"
#include "general/NumericalRecipes.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/ImageProcessing.hpp"
#include "process_main/ExStar.hpp"
#include "process_main/LinearSolve.hpp"
#include <mpi.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <array>
#include <cstddef>
#include <complex>
#include <memory>
#include <limits>
#include <new>
#include <stdexcept>

namespace PSFModel {

    using Internal::ChipPSFState;
    using Internal::CandidatePowerStatus;
    using Internal::ExposurePSFState;

    // Forward declarations of local helper functions
    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void starSelection(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void applyInitialFits(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);

    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio);
    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio);
    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio);
    void getPSFFWHM(
        const std::vector<float>& power,
        double& FWHM,
        int& star_area);

    // ==========================================
    // Function: Validate one PSF fitting sample
    // Method: Require finite position, shape diagnostics, and every power-spectrum pixel before retaining the star.
    // ==========================================
    static bool isFinitePSFStar(const std::vector<float>& star, std::size_t offset, int pixel_count,
                                double x, double y, double size, double e1, double e2) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(size) ||
            !std::isfinite(e1) || !std::isfinite(e2) ||
            offset + static_cast<std::size_t>(pixel_count) > star.size()) {
            return false;
        }
        for (int idx = 0; idx < pixel_count; ++idx) {
            if (!std::isfinite(star[offset + idx])) {
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Run Stage-5 Lite PSF modeling for one exposure
    // Method: Execute dynamic candidate loading, selection, diagnostics, and
    //         only the retained local-polynomial fitting path.
    // ==========================================
    void procPSF(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
        std::vector<std::string> imageFiles;
        std::string dirOutput;
        UniversalUtils::getImageList(expo_file_path, imageFiles, dirOutput);

        int nchip = static_cast<int>(imageFiles.size());

        auto state_ptr = std::make_unique<ExposurePSFState>(nchip);
        ExposurePSFState& state = *state_ptr;

        int nc = 0;
        std::vector<std::array<double, 4>> p_chip(
            static_cast<std::size_t>(nchip), {0.0, 0.0, 0.0, 0.0});

        readInCandidates(nchip, imageFiles, dirOutput, nc, p_chip, state);

        starSelection(nchip, imageFiles, dirOutput, state);

        applyInitialFits(nchip, imageFiles, dirOutput, state);

        plotStarExpo(nchip, imageFiles, dirOutput, state);

        plotStars(nchip, imageFiles, dirOutput, nc, p_chip, state);

        makePSFLocalFit(nchip, imageFiles, dirOutput, state);
    }

    // Local Helper Routines

    // ==========================================
    // Function: Load one exposure's PSF candidates and power stamps
    // Method: Gate each chip with the shared norm sentinel, reject numerically
    //         invalid corrected spectra, and compare only active candidates.
    // ==========================================
    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        const int ns = LensingConfig::ns;
        const int len_s = LensingConfig::len_s;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string headname = dirOutput + "/astrometry/Head/" + prefix_e + ".head";

        nc = 0;

        for (int k = 0; k < nchip; ++k) {
            ChipPSFState& chip = state.chips[k];
            chip.stars.clear();
            chip.stars.reserve(LensingConfig::nstar_max);
            chip.selection.clear();
            chip.fit.clear();

            const Universalblock::NormStatus norm_status =
                Universalblock::checkNorm(imageFiles[k], dirOutput);
            if (norm_status == Universalblock::NormStatus::Invalid) {
                continue;
            }
            if (norm_status != Universalblock::NormStatus::Valid) {
                MPIFailure::abortWorld(
                    "validate PSF chip norm",
                    Universalblock::normErrorDetail(
                        norm_status, imageFiles[k], dirOutput));
            }

            double cRPIX[2] = {0.0, 0.0};
            double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
            double cRVAL[2] = {0.0, 0.0};
            double PU[2][LensingConfig::npd];
            int ierror = 0;

            Astrometry::readAstrometryPara(headname, k + 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd, ierror);

            if (ierror == 1) continue;

            nc++;
            double x = 1.0;
            double y = 1.0;
            double xx = 0.0, yy = 0.0;
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][0] = xx;
            p_chip[nc - 1][1] = yy;

            x = 2046.0;
            y = 4094.0;
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][2] = xx;
            p_chip[nc - 1][3] = yy;

            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/dat_StarCanInfo", prefix, "_star_can_info.dat");

            std::ifstream infile(filepath);
            if (!infile.is_open()) {
                MPIFailure::abortWorld("read PSF star-candidate info", filepath);
            }

            std::string header;
            if (!std::getline(infile, header)) {
                MPIFailure::abortWorld(
                    "read PSF star-candidate header", filepath);
            }

            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                float aa[4] = {0.0f, 0.0f, 0.0f, 0.0f};
                if (iss >> aa[0] >> aa[1] >> aa[2] >> aa[3]) {
                    ChipPSFState::StarRow row{};
                    row[0] = aa[0];
                    row[1] = aa[1];
                    row[2] = aa[2];
                    row[3] = aa[3];
                    chip.stars.push_back(row);
                }
            }
            infile.close();

            const int nstar = state.getNStar(k);
            chip.selection.assign(
                static_cast<std::size_t>(nstar), Internal::StarSelectionState{});
            std::cout << "PSF candidates: chip=" << (k + 1)
                      << " count=" << nstar
                      << " selection_entries=" << chip.selection.size() << std::endl;

            if (nstar > 0) {
                int nn1 = ns * len_s;
                int nn2 = ns * ((nstar / len_s) + 1);
                std::string stampPath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                std::vector<float> star;
                if (!FitsIO::readStamps(nstar, 1, nstar, ns, ns, star, nn1, nn2, stampPath)) {
                    MPIFailure::abortWorld(
                        "read PSF star-candidate power", stampPath);
                }

                for (int i = 0; i < nstar; ++i) {
                    state.getStarPara(k, i, 4) = 1.0;

                    std::vector<float> source_p(ns * ns);
                    for (int v = 0; v < ns; ++v) {
                        for (int u = 0; u < ns; ++u) {
                            source_p[v * ns + u] = star[static_cast<size_t>(i) * ns * ns + v * ns + u];
                        }
                    }

                    double sum_power = 0.0;
                    double chi_window_sum = 0.0;
                    if (Internal::assessF77CandidatePower(
                            ns, ns, source_p, sum_power, chi_window_sum)
                        != CandidatePowerStatus::Accepted) {
                        state.getStarPara(k, i, 4) = -1.0;
                        continue;
                    }

                    x = state.getStarPara(k, i, 1);
                    y = state.getStarPara(k, i, 2);
                    Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
                    state.getStarPara(k, i, 5) = xx;
                    state.getStarPara(k, i, 6) = yy;

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, source_p, ee, size, 0.02f);
                    state.getStarPara(k, i, 7) = size;
                    state.getStarPara(k, i, 8) = ee[0];
                    state.getStarPara(k, i, 9) = ee[1];

                    double FWHM = 0.0;
                    int star_area = 0;
                    getPSFFWHM(source_p, FWHM, star_area);
                    if (!Internal::candidateDiagnosticsAreFinite(
                            size, ee[0], ee[1])) {
                        state.getStarPara(k, i, 4) = -1.0;
                        continue;
                    }
                    state.getStarPara(k, i, 10) = FWHM;
                    state.getStarPara(k, i, 11) = 1.0 / sum_power;
                    state.getStarPara(
                        k, i, ChipPSFState::star_area_index) =
                            static_cast<double>(star_area);

                    Internal::StarSelectionState& selection = chip.selection[i];
                    selection.full_power_sum = sum_power;
                    selection.chi_window_sum = chi_window_sum;
                    const Internal::PSFChiWindow window =
                        Internal::getPSFChiWindow(ns);
                    selection.chi_window.reserve(
                        static_cast<std::size_t>(window.pixelCount()));
                    for (int row = window.first; row <= window.last; ++row) {
                        for (int column = window.first; column <= window.last; ++column) {
                            selection.chi_window.push_back(source_p[row * ns + column]);
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // Function: Reject every PSF candidate in an exposure
    // Method: Clear F77 selection flags and release cached comparison windows.
    // ==========================================
    static void rejectExposureCandidates(ExposurePSFState& state) {
        for (ChipPSFState& chip : state.chips) {
            chip.fit.clear();
            for (std::size_t index = 0; index < chip.stars.size(); ++index) {
                chip.stars[index][4] = -1.0;
                if (index >= chip.selection.size()) continue;
                Internal::StarSelectionState& selection = chip.selection[index];
                selection.selected_group = false;
                selection.selected_fit = false;
                selection.min_chi = 1000.0f;
                selection.leverage = 0.0;
                std::vector<float>().swap(selection.chi_window);
            }
        }
    }

    // ==========================================
    // Function: Select PSF stars with the historical F77 algorithm
    // Method: Normalize every numerically safe candidate, form one exposure-wide
    //         chi threshold, and retain only each chip's largest valid group.
    // ==========================================
    static void starSelectionF77(int nchip, ExposurePSFState& state) {
        std::vector<std::vector<Internal::F77PSFCandidateView>> candidates(
            static_cast<std::size_t>(nchip));
        int safe_candidate_count = 0;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            int invalid_candidate_count = 0;
            candidates[chip_index].reserve(chip.selection.size());
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                Internal::StarSelectionState& selection =
                    chip.selection[star_index];
                selection.selected_group = false;
                selection.selected_fit = false;
                selection.min_chi = 1000.0f;
                selection.leverage = 0.0;
                Internal::F77PSFCandidateView candidate;
                if (!Internal::prepareF77PSFCandidate(
                        star_index,
                        state.getStarPara(chip_index, star_index, 4),
                        selection.full_power_sum,
                        state.getStarPara(chip_index, star_index, 7),
                        selection.chi_window,
                        candidate)) {
                    state.getStarPara(chip_index, star_index, 4) = -1.0;
                    std::vector<float>().swap(selection.chi_window);
                    invalid_candidate_count++;
                    continue;
                }
                candidates[chip_index].push_back(candidate);
                safe_candidate_count++;
            }
            std::cout << "PSF_F77_CHIP chip=" << (chip_index + 1)
                      << " invalid_candidates=" << invalid_candidate_count
                      << " safe_candidates=" << candidates[chip_index].size()
                      << '\n';
        }

        if (!Internal::hasMinimumF77PSFCandidates(
                static_cast<std::size_t>(safe_candidate_count),
                LensingConfig::nstar_min)) {
            std::cout << "PSF_F77_EXPOSURE candidates=" << safe_candidate_count
                      << " decision=REJECT_MINIMUM\n";
            rejectExposureCandidates(state);
            return;
        }

        Internal::F77PSFPairStatistics statistics;
        if (!Internal::computeF77PSFPairStatistics(candidates, statistics)
            || statistics.threshold_pair_chi.size() <= 4U) {
            std::cout << "PSF_F77_EXPOSURE decision=REJECT_INVALID_PAIRS\n";
            rejectExposureCandidates(state);
            return;
        }
        float peak = 0.0f;
        float width = 0.0f;
        NumericalRecipes::getPeakWidthLowSide(
            statistics.threshold_pair_chi, peak, width);
        const float chi_threshold = peak + 4.0f * width;
        if (!std::isfinite(peak) || !std::isfinite(width)
            || !std::isfinite(chi_threshold)) {
            std::cout << "PSF_F77_EXPOSURE decision=REJECT_INVALID_THRESHOLD\n";
            rejectExposureCandidates(state);
            return;
        }

        std::vector<std::vector<int>> selected;
        if (!Internal::selectF77PSFLargestGroups(
                candidates, statistics, chi_threshold,
                LensingConfig::nstar_min_local, selected)) {
            std::cout << "PSF_F77_EXPOSURE decision=REJECT_INVALID_GRAPH\n";
            rejectExposureCandidates(state);
            return;
        }
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            std::vector<bool> keep(
                static_cast<std::size_t>(state.getNStar(chip_index)), false);
            for (int star_index : selected[chip_index]) keep[star_index] = true;
            for (std::size_t compact_index = 0;
                 compact_index < candidates[chip_index].size();
                 ++compact_index) {
                const int original_index =
                    candidates[chip_index][compact_index].star_index;
                chip.selection[original_index].min_chi =
                    statistics.min_chi[chip_index][compact_index];
            }
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                Internal::StarSelectionState& selection =
                    chip.selection[star_index];
                selection.selected_group = keep[star_index];
                state.getStarPara(chip_index, star_index, 4) =
                    keep[star_index] ? 1.0 : -1.0;
                if (!keep[star_index]) {
                    std::vector<float>().swap(selection.chi_window);
                }
            }
            std::cout << "PSF_F77_CHIP chip=" << (chip_index + 1)
                      << " selected=" << selected[chip_index].size()
                      << " threshold=" << chi_threshold << '\n';
        }
    }

    // ==========================================
    // Function: Select PSF stars with the frozen Lite policy
    // Method: Run the sole F77-compatible exposure grouping implementation.
    // ==========================================
    void starSelection(
        int nchip,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        ExposurePSFState& state) {
        (void)imageFiles;
        (void)dirOutput;
        starSelectionF77(nchip, state);
    }

    // ==========================================
    // Structure: Hold one chip's ordered polynomial fitting samples
    // Method: Preserve original candidate indices alongside positions and
    //         contiguous row-major power stamps.
    // ==========================================
    struct ChipFitSamples {
        std::vector<int> star_indices;
        std::vector<std::array<double, 2>> positions;
        std::vector<float> power;
    };

    // ==========================================
    // Function: Read all candidate power stamps for one chip
    // Method: Use the live candidate count and existing sharded Stage-4 product.
    // ==========================================
    static std::vector<float> readChipCandidatePower(
        int chip_index,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        const ExposurePSFState& state) {
        const int nstar = state.getNStar(chip_index);
        if (nstar <= 0) return {};
        const int ns = LensingConfig::ns;
        const int nn1 = ns * LensingConfig::len_s;
        const int nn2 = ns * ((nstar / LensingConfig::len_s) + 1);
        const std::string prefix = UniversalUtils::getPrefix(imageFiles[chip_index]);
        const std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/fits_StarCanP", prefix,
            "_star_can_power.fits");
        std::vector<float> power;
        if (!FitsIO::readStamps(
                nstar, 1, nstar, ns, ns, power, nn1, nn2, filename)) {
            MPIFailure::abortWorld("read initial-fit PSF star power", filename);
        }
        return power;
    }

    // ==========================================
    // Function: Gather an ordered subset of one chip's fitting samples
    // Method: Validate each requested original candidate and copy its position
    //         and full power stamp into a compact fitting buffer.
    // ==========================================
    static bool buildChipFitSamples(
        int chip_index,
        const std::vector<int>& selected_indices,
        const std::vector<float>& all_power,
        const ExposurePSFState& state,
        ChipFitSamples& samples) {
        samples = {};
        const int ns = LensingConfig::ns;
        const int pixel_count = ns * ns;
        samples.star_indices.reserve(selected_indices.size());
        samples.positions.reserve(selected_indices.size());
        samples.power.reserve(selected_indices.size()
                              * static_cast<std::size_t>(pixel_count));
        for (int star_index : selected_indices) {
            if (star_index < 0 || star_index >= state.getNStar(chip_index)) return false;
            const double x = state.getStarPara(chip_index, star_index, 1);
            const double y = state.getStarPara(chip_index, star_index, 2);
            const double size = state.getStarPara(chip_index, star_index, 7);
            const double e1 = state.getStarPara(chip_index, star_index, 8);
            const double e2 = state.getStarPara(chip_index, star_index, 9);
            const std::size_t offset =
                static_cast<std::size_t>(star_index) * pixel_count;
            if (!isFinitePSFStar(
                    all_power, offset, pixel_count, x, y, size, e1, e2)) {
                return false;
            }
            samples.star_indices.push_back(star_index);
            samples.positions.push_back({x, y});
            samples.power.insert(
                samples.power.end(),
                all_power.begin() + static_cast<std::ptrdiff_t>(offset),
                all_power.begin() + static_cast<std::ptrdiff_t>(offset + pixel_count));
        }
        return true;
    }

    // ==========================================
    // Function: Fit one ordered chip sample set and compute leverage
    // Method: Use Lite's sole normalized local-polynomial design and share the
    //         fitted coefficients and hat diagonals with final output.
    // ==========================================
    static LinearSolve::SolveStatus fitChipSamples(
        const ChipFitSamples& samples,
        std::vector<double>& coefficients,
        std::vector<double>& leverage,
        LinearSolve::SolveDiagnostics& diagnostics) {
        const int sample_count = static_cast<int>(samples.star_indices.size());
        return itpNormPSF(
            sample_count, samples.power, samples.positions,
            LensingConfig::ns, LensingConfig::npl,
            LensingConfig::chipnx, LensingConfig::chipny,
            coefficients, &diagnostics, &leverage);
    }

    // ==========================================
    // Function: Handle a failed initial PSF fit
    // Method: Preserve F77 group flags and invalidate only the unusable fit cache.
    // ==========================================
    static void handleInitialFitFailure(
        int chip_index,
        ExposurePSFState& state) {
        state.chips[chip_index].fit.clear();
    }

    // ==========================================
    // Function: Build the initial selected-star fit cache for one chip
    // Method: Fit the F77 grouping survivors once and cache coefficients and
    //         leverage without applying any rejection or refit policy.
    // ==========================================
    static bool buildInitialFitCache(
        int chip_index,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        ExposurePSFState& state) {
        ChipPSFState& chip = state.chips[chip_index];
        chip.fit.clear();
        std::vector<int> selected_indices;
        for (int star_index = 0;
             star_index < state.getNStar(chip_index); ++star_index) {
            Internal::StarSelectionState& selection =
                chip.selection[star_index];
            selection.selected_fit = selection.selected_group;
            selection.leverage = 0.0;
            if (selection.selected_group) selected_indices.push_back(star_index);
        }
        if (static_cast<int>(selected_indices.size())
            < LensingConfig::nstar_min_local) {
            handleInitialFitFailure(chip_index, state);
            return false;
        }

        const std::vector<float> all_power = readChipCandidatePower(
            chip_index, imageFiles, dirOutput, state);
        ChipFitSamples samples;
        if (!buildChipFitSamples(
                chip_index, selected_indices, all_power, state, samples)) {
            LinearSolve::reportFailure(
                "PSFModel::buildInitialFitCache",
                LinearSolve::SolveStatus::FailedSolver,
                "chip=" + std::to_string(chip_index + 1)
                    + " reason=INVALID_SELECTED_SAMPLE"
                      " action=KEEP_F77_SELECTION_UNFITTED");
            handleInitialFitFailure(chip_index, state);
            return false;
        }

        LinearSolve::SolveDiagnostics diagnostics;
        std::vector<double> coefficients;
        std::vector<double> leverage;
        const LinearSolve::SolveStatus status = fitChipSamples(
            samples, coefficients, leverage, diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            LinearSolve::reportFailure(
                "PSFModel::buildInitialFitCache", status,
                "chip=" + std::to_string(chip_index + 1) + " "
                    + LinearSolve::diagnosticsContext(diagnostics)
                    + " action=KEEP_F77_SELECTION_UNFITTED");
            handleInitialFitFailure(chip_index, state);
            return false;
        }

        chip.fit.valid = true;
        chip.fit.initial_star_count =
            static_cast<int>(samples.star_indices.size());
        chip.fit.star_indices = samples.star_indices;
        chip.fit.coefficients = std::move(coefficients);
        chip.fit.leverage = std::move(leverage);
        for (std::size_t local_index = 0;
             local_index < chip.fit.star_indices.size(); ++local_index) {
            const int star_index = chip.fit.star_indices[local_index];
            chip.selection[star_index].selected_fit = true;
            chip.selection[star_index].leverage =
                chip.fit.leverage[local_index];
        }
        return true;
    }

    // ==========================================
    // Function: Build all frozen Lite Stage-5 initial fit caches
    // Method: Fit each chip once after F77 grouping and skip rejection/refitting.
    // ==========================================
    void applyInitialFits(
        int nchip,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        ExposurePSFState& state) {
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const bool fitted = buildInitialFitCache(
                chip_index, imageFiles, dirOutput, state);
            std::cout << "PSF_F77_INITIAL_FIT chip=" << (chip_index + 1)
                      << " selected="
                      << state.chips[chip_index].fit.initial_star_count
                      << " decision="
                      << (fitted ? "KEEP_INITIAL_FIT"
                                 : "KEEP_F77_SELECTION_UNFITTED")
                      << '\n';
        }
    }

    // ==========================================
    // Function: Assemble exposure-wide selected-star power stamps
    // Method: Append each live chip buffer and retain every selected star while
    //         preserving the legacy two-dimensional FITS mosaic layout.
    // ==========================================
    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        const int ns = LensingConfig::ns;
        const int len_s = LensingConfig::len_s;
        std::vector<int> opt;
        std::vector<float> star_test;
        int ntot = 0;
        int start = 0;

        for (int ichip = 0; ichip < nchip; ++ichip) {
            const int nstar = state.getNStar(ichip);
            if (nstar == 0) continue;
            int nn1 = ns * len_s;
            int nn2 = ns * ((nstar / len_s) + 1);

            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", prefix,
                "_star_can_power.fits");
            std::vector<float> star;
            if (!FitsIO::readStamps(
                    nstar, 1, nstar, ns, ns, star, nn1, nn2, filepath)) {
                MPIFailure::abortWorld(
                    "read exposure PSF star power", filepath);
            }
            star_test.insert(star_test.end(), star.begin(), star.end());
            opt.resize(static_cast<std::size_t>(start + nstar), 0);

            for (int i = 0; i < nstar; ++i) {
                if (state.getStarPara(ichip, i, 4) <= 0.0) continue;
                ntot++;
                opt[start + i] = 1;
            }
            start += nstar;
        }

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string out_filename = dirOutput + "/stamps/fits_StarP/"
            + prefix_e + "_star_power_expo.fits";
        if (ntot > 0) {
            int len_sam = LensingConfig::len_sam;
            int nn1 = ns * len_sam;
            int nn2 = ns * ((ntot / len_sam) + 1);
            FitsIO::writeStamps2(
                start, start, ns, ns, star_test, opt, 1, nn1, nn2,
                out_filename);
        }
    }

    // ==========================================
    // Function: Publish exposure-level star diagnostics and PSF visualization
    // Method: Aggregate selected stars and route text/FITS products through
    //         checked main-process writers.
    // ==========================================
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        (void)nc;
        int nm = 1000;
        int nstar_min_local = LensingConfig::nstar_min_local;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        // prefix_dir inlined: per-type stamps/ subdirs (reorganized layout)

        std::string info_filename = dirOutput + "/stamps/dat_StarInfo/" + prefix_e + "_star_info_expo.dat";
        MainIO::OutputFile outfile(info_filename);
        if (!outfile.is_open()) {
            std::cerr << "plotStars: Error opening " << info_filename << std::endl;
            return;
        }
        outfile << std::setprecision(10);

        outfile << "# ichip nstar FWHM e1 e2 chi_d\n";

        std::vector<std::array<double, 5>> sk;
        std::size_t total_candidates = 0;
        for (int k = 0; k < nchip; ++k) {
            total_candidates += static_cast<std::size_t>(state.getNStar(k));
        }
        sk.reserve(total_candidates);

        for (int k = 0; k < nchip; ++k) {
            double FWHM_ave = 0.0;
            double e1_ave = 0.0;
            double e2_ave = 0.0;
            double chi_d_ave = 0.0;
            int nums = 0;
            int prev_idx = -1;

            for (int i = 0; i < state.getNStar(k); ++i) {
                if (state.getStarPara(k, i, 4) <= 0.0) continue;
                nums++;
                sk.push_back({
                    state.getStarPara(k, i, 5),
                    state.getStarPara(k, i, 6),
                    state.getStarPara(k, i, 7),
                    state.getStarPara(k, i, 8),
                    state.getStarPara(k, i, 9)
                });

                FWHM_ave += state.getStarPara(k, i, 10);
                e1_ave += state.getStarPara(k, i, 8);
                e2_ave += state.getStarPara(k, i, 9);
                if (nums >= 2 && prev_idx != -1) {
                    chi_d_ave += Internal::normalizedChiDistance(
                        state.chips[k].selection[i].chi_window,
                        state.chips[k].selection[prev_idx].chi_window);
                }
                prev_idx = i;
            }

            if (nums >= nstar_min_local) {
                FWHM_ave /= nums;
                e1_ave /= nums;
                e2_ave /= nums;
                chi_d_ave /= (nums - 1.0);
                outfile << (k + 1) << " " << nums << " "
                        << std::scientific << std::setprecision(10)
                        << FWHM_ave << " " << e1_ave << " " << e2_ave << " " << chi_d_ave << "\n";
            } else {
                outfile << (k + 1) << " 0 -99.0 -99.0 -99.0 -99.0\n";
            }
        }

        // std::cout << prefix_dir << " total no. of stars: " << sk.size() << std::endl;
        outfile.close();

        std::vector<float> PSFmap;
        ImageProcessing::drawShearExpo(nm, PSFmap, p_chip, sk, 200.0, 1.0);

        std::string fits_filename = dirOutput + "/stamps/fits_PsfSrc/" + prefix_e + "_PSF_source.fits";
        FitsIO::writeImage(fits_filename, nm, nm, PSFmap);
    }

    // ==========================================
    // Function: Fit and serialize local PSF models.
    // Method: Preserve F77 model layout with 17-digit serialization while
    //         reporting the ordinary full-fit model for each fitted star.
    // ==========================================
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int len_s = LensingConfig::len_s;
        int npl = LensingConfig::npl;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string comp_filename = dirOutput + "/stamps/dat_StarComp/" + prefix_e + "_star_comp_expo.dat";
        MainIO::OutputFile file90(comp_filename);
        if (!file90.is_open()) {
            std::cerr << "makePSFLocalFit: Error opening " << comp_filename << std::endl;
            return;
        }
        file90 << std::setprecision(17);

        for (int k = 0; k < nchip; ++k) {
            int nums = 0;
            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);

            std::vector<float> star;
            if (state.getNStar(k) > 0) {
                int nn1 = ns * len_s;
                int nn2 = ns * ((state.getNStar(k) / len_s) + 1);
                std::string filepath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                if (!FitsIO::readStamps(
                        state.getNStar(k), 1, state.getNStar(k), ns, ns,
                        star, nn1, nn2, filepath)) {
                    MPIFailure::abortWorld(
                        "read local-fit PSF star power", filepath);
                }
            }

            std::string coe_filename = OutputLayout::chipPath(
                dirOutput, "stamps/dat_PsfFit", prefix, "_PSF_coe_local.dat");
            MainIO::OutputFile file10(coe_filename);
            if (!file10.is_open()) {
                std::cerr << "makePSFLocalFit: Error opening " << coe_filename << std::endl;
                continue;
            }
            file10 << std::setprecision(17);

            std::vector<std::array<double, 2>> posi;
            std::vector<std::array<double, 3>> sshape;
            int removed_non_finite = 0;

            const Internal::ChipPSFFitState& cached_fit = state.chips[k].fit;
            for (int i : cached_fit.star_indices) {
                double px = state.getStarPara(k, i, 1);
                double py = state.getStarPara(k, i, 2);
                double shape_size = state.getStarPara(k, i, 7);
                double shape_e1 = state.getStarPara(k, i, 8);
                double shape_e2 = state.getStarPara(k, i, 9);
                std::size_t star_offset = static_cast<std::size_t>(i) * ns * ns;
                if (!isFinitePSFStar(star, star_offset, ns * ns, px, py,
                                     shape_size, shape_e1, shape_e2)) {
                    removed_non_finite++;
                    continue;
                }
                nums++;
                posi.push_back({px, py});
                sshape.push_back({shape_size, shape_e1, shape_e2});
            }
            std::vector<double> PSF_coe_l = cached_fit.coefficients;
            LinearSolve::SolveDiagnostics fit_diagnostics;
            LinearSolve::SolveStatus fit_status = cached_fit.valid
                ? LinearSolve::SolveStatus::Normal
                : LinearSolve::SolveStatus::FailedRankDeficient;

            if (nums >= LensingConfig::nstar_min_local &&
                fit_status == LinearSolve::SolveStatus::Normal &&
                cached_fit.star_indices.size()
                    == static_cast<std::size_t>(nums)) {

                file90 << (k + 1) << " " << nums << " 1\n";

                std::vector<float> poly_cochi2(nums);
                float poly_ave = 0.0f, poly_std = 0.0f;

                for (int i = 0; i < nums; ++i) {
                    double xx = 2.0 * (posi[i][0] / static_cast<double>(LensingConfig::chipnx)) - 1.0;
                    double yy = 2.0 * (posi[i][1] / static_cast<double>(LensingConfig::chipny)) - 1.0;
                    std::vector<float> model, model0;
                    getPSFModel(ns, npl, PSF_coe_l, xx, yy, model, model0);
                    ExStar::anaChi2Simple(ns, model.data(), model0.data(), poly_cochi2[i]);

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, model, ee, size, 0.02f);

                    double msshape_size = size;
                    double msshape_e1 = ee[0];
                    double msshape_e2 = ee[1];

                    float px = static_cast<float>(posi[i][0]);
                    float py = static_cast<float>(posi[i][1]);

                    file90 << px << " " << py << " "
                           << sshape[i][0] << " " << sshape[i][1] << " " << sshape[i][2] << " "
                           << msshape_size << " " << msshape_e1 << " " << msshape_e2 << "\n";
                }

                ExStar::getArrayAveStd(poly_cochi2, poly_ave, poly_std);
                file10 << nums << " 1 " << poly_ave << " " << poly_std << "\n";
                for (int i = 0; i < ns; ++i) {
                    for (int j = 0; j < ns; ++j) {
                        for (int u = 0; u < npl + 1; ++u) {
                            file10 << PSF_coe_l[(j * ns + i) * (npl + 1) + u] << (u == npl ? "" : " ");
                        }
                        file10 << "\n";
                    }
                }
            } else {
                if (nums < LensingConfig::nstar_min_local) {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", LinearSolve::SolveStatus::FailedRankDeficient,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " valid_samples=" + std::to_string(nums) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " required=" + std::to_string(LensingConfig::nstar_min_local) +
                            " action=MARK_CHIP_INVALID");
                } else {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", fit_status,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " " + LinearSolve::diagnosticsContext(fit_diagnostics) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " action=MARK_CHIP_INVALID");
                }
                file10 << nums << " -1 -1 -1\n";
                file90 << (k + 1) << " " << nums << " -1\n";
            }
            file10.close();
        }
        file90.close();
    }

    // Mathematical Interpolation Helpers

    // ==========================================
    // Function: Fit all PSF-frequency pixels with a shared spatial design
    // Method: Factor the constant and polynomial designs once, then reuse both QR objects for every RHS pixel.
    // Note:   Lite retains only the normalized local-polynomial design, so the
    //         chip-coordinate mapping to [-1,1] is unconditional.
    // ==========================================
    static LinearSolve::SolveStatus fitPSFCoefficients(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny,
        std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage) {
        const int pixel_count = ns * ns;
        PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);

        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = nsam;
        diag.cols = npp;
        diag.required_rank = npp;

        if (nsam < npp || ns <= 0 || npp <= 0 ||
            static_cast<int>(posi.size()) < nsam ||
            image.size() < static_cast<size_t>(nsam) * pixel_count ||
            nx <= 0 || ny <= 0) {
            diag.rank = std::min(nsam, npp);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd polynomial_design(nsam, npp);
        Eigen::MatrixXd constant_design = Eigen::MatrixXd::Ones(nsam, 1);
        for (int sample = 0; sample < nsam; ++sample) {
            double xx = 2.0 * (posi[sample][0] / static_cast<double>(nx)) - 1.0;
            double yy = 2.0 * (posi[sample][1] / static_cast<double>(ny)) - 1.0;
            if (!std::isfinite(xx) || !std::isfinite(yy)) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
            for (int term = 0; term < npp; ++term) {
                polynomial_design(sample, term) = UniversalUtils::fitFunc2(xx, yy, term);
            }
        }
        if (!polynomial_design.allFinite()) {
            return LinearSolve::SolveStatus::FailedSolver;
        }
        for (std::size_t idx = 0; idx < static_cast<std::size_t>(nsam) * pixel_count; ++idx) {
            if (!std::isfinite(image[idx])) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
        }

        LinearSolve::LeastSquaresQR constant_solver;
        LinearSolve::SolveDiagnostics constant_diagnostics;
        LinearSolve::SolveStatus status = constant_solver.factorize(constant_design, constant_diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            diag = constant_diagnostics;
            return status;
        }

        LinearSolve::LeastSquaresQR polynomial_solver;
        status = polynomial_solver.factorize(polynomial_design, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        if (leverage != nullptr) {
            Eigen::MatrixXd covariance_base;
            status = polynomial_solver.unscaledCovariance(covariance_base);
            if (status != LinearSolve::SolveStatus::Normal) return status;
            leverage->assign(static_cast<std::size_t>(nsam), 0.0);
            for (int sample = 0; sample < nsam; ++sample) {
                const double value = (
                    polynomial_design.row(sample) * covariance_base
                    * polynomial_design.row(sample).transpose())(0, 0);
                double loo_residual = 0.0;
                double loo_model = 0.0;
                if (!Internal::computeAnalyticLOO(
                        0.0, 0.0, value,
                        LensingConfig::psf_loo_min_denom,
                        loo_residual, loo_model)) {
                    leverage->clear();
                    return LinearSolve::SolveStatus::FailedIllConditioned;
                }
                (*leverage)[sample] = std::max(0.0, value);
            }
        }

        Eigen::VectorXd rhs(nsam);
        Eigen::VectorXd constant_solution;
        Eigen::VectorXd polynomial_solution;
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                const int pixel = j * ns + i;
                for (int sample = 0; sample < nsam; ++sample) {
                    rhs(sample) = image[static_cast<size_t>(sample) * pixel_count + pixel];
                }

                status = constant_solver.solve(rhs, constant_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }
                status = polynomial_solver.solve(rhs, polynomial_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }

                for (int term = 0; term < npp; ++term) {
                    PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + term] = polynomial_solution(term);
                }
                PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + npp] = constant_solution(0);
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // ==========================================
    // Function: Fit normalized local PSF coefficients
    // Method: Normalize chip coordinates to [-1,1] and delegate to the shared-design batch fitter.
    // ==========================================
    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage) {
        return fitPSFCoefficients(
            nsam, image, posi, ns, npp, nx, ny,
            PSF_coe, diagnostics, leverage);
    }

    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0) {
        modelp.assign(static_cast<size_t>(ns) * ns, 0.0f);
        model0.assign(static_cast<size_t>(ns) * ns, 0.0f);
        std::vector<double> coep(npp);
        std::vector<double> coe0(1);
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                for (int k = 0; k < npp; ++k) {
                    coep[k] = PSF_coe[(j * ns + i) * (npp + 1) + k];
                }
                coe0[0] = PSF_coe[(j * ns + i) * (npp + 1) + npp];
                model0[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, 1, coe0));
                modelp[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, npp, coep));
                if (std::isnan(modelp[j * ns + i])) {
                    modelp[0] = modelp[j * ns + i];
                    return;
                }
            }
        }
    }

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;
                        }
                    }
                }
            }
            area0 = tempi;
        }

        area = (area_cnt - 1) / 2;
    }

    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        e[0] = 0.0;
        e[1] = 0.0;
        double norm = 0.0;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;

                            double kx = u - cx;
                            double ky = v - cy;
                            double p_val = power[idx];
                            e[0] += p_val * (kx * kx - ky * ky);
                            e[1] += p_val * 2.0 * kx * ky;
                            norm += p_val * (kx * kx + ky * ky);
                        }
                    }
                }
            }
            area0 = tempi;
        }

        e[0] /= norm;
        e[1] /= norm;
    }

    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio) {
        int area = 0;
        getPowerArea(nx, ny, power, area, thresh_ratio);
        size = area;
        getPowerE(nx, ny, power, e, thresh_ratio);
    }

    // ==========================================
    // Function: Measure historical PSF FWHM and exact exp(-1) star area
    // Method: Count integer threshold pixels once, then reuse the unchanged
    //         area-minus-1e-5 conversion for the reported FWHM.
    // ==========================================
    void getPSFFWHM(
        const std::vector<float>& power,
        double& FWHM,
        int& star_area) {
        star_area = Internal::countPSFStarArea(power, LensingConfig::ns);
        FWHM = Internal::fwhmFromStarArea(
            static_cast<double>(star_area),
            LensingConfig::ns,
            LensingConfig::pixel_size);
    }
}
