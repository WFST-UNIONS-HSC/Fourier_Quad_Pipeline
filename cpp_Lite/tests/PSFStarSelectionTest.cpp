#include "process_main/PSFStarSelection.hpp"
#include "process_main/PSFModelState.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace {

using namespace PSFModel::Internal;

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF star-selection test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void testChiWindowAndDistance() {
    const PSFChiWindow window = getPSFChiWindow(64);
    require(window.first == 15 && window.last == 47
                && window.pixelCount() == 1089,
            "64x64 chi window must remain inclusive 15..47");
    std::vector<float> first(9, 0.1f);
    std::vector<float> second = first;
    require(normalizedChiDistance(first, second) == 0.0f,
            "identical windows must have zero chi distance");
    second[4] += 0.05f;
    require(normalizedChiDistance(first, second) > 0.0f,
            "different windows must have positive chi distance");
}

void testStarAreaMeasurementAndStorage() {
    std::vector<float> power(25, 0.0f);
    power[12] = 10.0f;
    power[6] = 4.0f;
    power[7] = 5.0f;
    power[11] = 6.0f;
    require(countPSFStarArea(power, 5) == 4,
            "star area must preserve the exp(-1) threshold count");
    const double fwhm = fwhmFromStarArea(4.0, 5, 0.2628);
    require(std::isfinite(fwhm) && fwhm > 0.0,
            "positive star area must yield finite FWHM");

    using ChipState = PSFModel::Internal::ChipPSFState;
    static_assert(ChipState::star_area_index == 12);
    static_assert(ChipState::star_area_index
        < static_cast<int>(std::tuple_size<ChipState::StarRow>::value));
    ChipState::StarRow row{};
    row[7] = 17.0;
    row[10] = fwhm;
    row[11] = 0.25;
    row[ChipState::star_area_index] = 4.0;
    require(row[7] == 17.0 && row[10] == fwhm && row[11] == 0.25
                && row[ChipState::star_area_index] == 4.0,
            "star-area storage must not alter existing fields");
}

void testF77CandidateSafetyPopulation() {
    auto countSafeCandidates = [](int raw_count) {
        std::vector<std::vector<float>> windows(
            static_cast<std::size_t>(raw_count), {0.8f, 0.2f});
        std::vector<F77PSFCandidateView> candidates;
        for (int index = 0; index < raw_count; ++index) {
            F77PSFCandidateView candidate;
            if (prepareF77PSFCandidate(
                    index, 1.0, index == 0 ? 0.0 : 1.0,
                    static_cast<double>(index + 1), windows[index], candidate)) {
                candidates.push_back(candidate);
            }
        }
        return candidates.size();
    };
    require(countSafeCandidates(21) == 20U,
            "one invalid candidate must not reject safe chip peers");
    require(!hasMinimumF77PSFCandidates(countSafeCandidates(192), 96)
                && hasMinimumF77PSFCandidates(countSafeCandidates(193), 96),
            "the F77 exposure minimum must count only safe candidates");

    std::vector<float> negative_sum_window = {-0.8f, -0.2f};
    F77PSFCandidateView candidate;
    require(prepareF77PSFCandidate(
                4, 1.0, -1.0, 2.0, negative_sum_window, candidate)
                && candidate.star_index == 4
                && negative_sum_window[0] > 0.0f,
            "finite negative sums must remain F77-compatible");

    std::vector<float> nonfinite_window = {
        0.8f, std::numeric_limits<float>::quiet_NaN()};
    require(!prepareF77PSFCandidate(
                5, 1.0, 1.0, 2.0, nonfinite_window, candidate),
            "non-finite windows must be rejected individually");
}

void testF77SelectionPolicy() {
    std::vector<std::vector<float>> rank_windows = {
        {0.80f, 0.20f}, {0.79f, 0.21f}, {0.78f, 0.22f},
        {0.77f, 0.23f}, {0.76f, 0.24f}, {0.75f, 0.25f}};
    std::vector<std::vector<F77PSFCandidateView>> rank_candidates(1);
    for (int index = 0; index < 6; ++index) {
        rank_candidates[0].push_back({
            index, static_cast<double>(index + 1), &rank_windows[index]});
    }
    F77PSFPairStatistics rank_statistics;
    require(computeF77PSFPairStatistics(rank_candidates, rank_statistics)
                && rank_statistics.size_threshold == 4.0
                && rank_statistics.threshold_pair_chi.size() == 3,
            "F77 size rank and both-large pair sampling must remain exact");

    std::vector<std::vector<float>> windows = {
        {0.80f, 0.20f}, {0.20f, 0.80f},
        {0.81f, 0.19f}, {0.82f, 0.18f}};
    const std::array<int, 4> original_indices = {0, 2, 5, 7};
    std::vector<std::vector<F77PSFCandidateView>> candidates(1);
    for (std::size_t compact = 0; compact < windows.size(); ++compact) {
        candidates[0].push_back({
            original_indices[compact], static_cast<double>(compact + 1),
            &windows[compact]});
    }
    F77PSFPairStatistics statistics;
    std::vector<std::vector<int>> selected;
    require(computeF77PSFPairStatistics(candidates, statistics)
                && selectF77PSFLargestGroups(
                    candidates, statistics, 0.05f, 2, selected),
            "finite F77 grouping must succeed");
    require(selected[0] == std::vector<int>({0, 5, 7}),
            "largest group must map compact entries to original indices");
    require(selectF77PSFLargestGroups(
                candidates, statistics, 0.05f, 4, selected)
                && selected[0].empty(),
            "a group below the local minimum must reject the chip");
}

void testAnalyticLOO() {
    double residual = 0.0;
    double model = 0.0;
    require(computeAnalyticLOO(
                10.0, 8.0, 0.25, 1.0e-6, residual, model)
                && std::abs(residual - 8.0 / 3.0) < 1.0e-12
                && std::abs(model - 22.0 / 3.0) < 1.0e-12,
            "analytic LOO must preserve the leverage correction");
    require(!computeAnalyticLOO(
                1.0, 0.0, 1.0 - 5.0e-7,
                1.0e-6, residual, model),
            "near-unit leverage must fail the denominator guard");
}

}  // namespace

int main() {
    testChiWindowAndDistance();
    testStarAreaMeasurementAndStorage();
    testF77CandidateSafetyPopulation();
    testF77SelectionPolicy();
    testAnalyticLOO();
    std::cout << "PSF star-selection tests passed\n";
    return EXIT_SUCCESS;
}
