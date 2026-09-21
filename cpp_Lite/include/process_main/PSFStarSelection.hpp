#ifndef PSF_STAR_SELECTION_HPP
#define PSF_STAR_SELECTION_HPP

#include <cstddef>
#include <vector>

namespace PSFModel {
namespace Internal {
// ==========================================
// Structure: Describe the inclusive central Fourier window used by PSF chi
// Method: Derive one shared bound pair from the square stamp side length.
// ==========================================
struct PSFChiWindow {
    int first = 0;
    int last = -1;

    // ==========================================
    // Function: Return the inclusive chi-window side length
    // Method: Convert valid bounds to a positive size and invalid bounds to zero.
    // ==========================================
    int side() const { return last >= first ? last - first + 1 : 0; }

    // ==========================================
    // Function: Return the square chi-window pixel count
    // Method: Square the derived inclusive side length.
    // ==========================================
    int pixelCount() const { return side() * side(); }
};

// ==========================================
// Function: Derive the shared PSF chi-window bounds
// Method: Preserve the legacy inclusive n/4-1 through 3n/4-1 definition.
// ==========================================
PSFChiWindow getPSFChiWindow(int n);
// ==========================================
// Function: Compute the normalized F77 PSF chi distance
// Method: Compare equal cached windows using the legacy signed-flux denominator.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second);

// ==========================================
// Function: Count exp(-1)-threshold pixels in one square Fourier-power stamp
// Method: Compare every finite stamp value with the central value times exp(-1).
// ==========================================
int countPSFStarArea(
    const std::vector<float>& power,
    int stamp_side);

// ==========================================
// Function: Convert an exp(-1) star area to the historical PSF FWHM
// Method: Preserve the legacy area-minus-1e-5 formula and explicit pixel scale.
// ==========================================
double fwhmFromStarArea(
    double star_area,
    int stamp_side,
    double pixel_size);
// Structure: View one numerically safe candidate in F77 selection
// Method: Preserve its original chip-local index while borrowing the normalized
//         central window and exposing its legacy size.
// ==========================================
struct F77PSFCandidateView {
    int star_index = -1;
    double size = 0.0;
    const std::vector<float>* chi_window = nullptr;
};

// ==========================================
// Function: Prepare one candidate for the compact F77 selection population
// Method: Apply candidate-level numerical gates, normalize its borrowed window,
//         and retain the original chip-local star index without a positivity gate.
// ==========================================
bool prepareF77PSFCandidate(
    int star_index,
    double selection_flag,
    double full_power_sum,
    double size,
    std::vector<float>& chi_window,
    F77PSFCandidateView& candidate);

// ==========================================
// Function: Test the historical F77 exposure candidate minimum
// Method: Compare only the compact numerically safe population with twice the
//         configured per-exposure minimum.
// ==========================================
bool hasMinimumF77PSFCandidates(
    std::size_t safe_candidate_count,
    int minimum_stars);

// ==========================================
// Structure: Store the single-pass F77 pair statistics
// Method: Preserve the exposure size rank, per-candidate 1000-sentinel minChi,
//         and the both-large-endpoint threshold sample.
// ==========================================
struct F77PSFPairStatistics {
    bool valid = false;
    double size_threshold = 0.0;
    std::vector<std::vector<float>> min_chi;
    std::vector<float> threshold_pair_chi;
};

// ==========================================
// Function: Compute the exact F77 exposure-size and same-chip pair statistics
// Method: Use the one-based floor(2N/3) rank, visit every unordered pair once,
//         and fail explicitly when a supposedly safe pair has non-finite chi.
// ==========================================
bool computeF77PSFPairStatistics(
    const std::vector<std::vector<F77PSFCandidateView>>& candidates_by_chip,
    F77PSFPairStatistics& statistics);

// ==========================================
// Function: Select only the largest F77 threshold component on every chip
// Method: Apply inclusive minChi/edge cuts, both local-minimum checks, and the
//         historical first-component tie behavior without secondary groups.
// ==========================================
bool selectF77PSFLargestGroups(
    const std::vector<std::vector<F77PSFCandidateView>>& candidates_by_chip,
    const F77PSFPairStatistics& statistics,
    float chi_threshold,
    int minimum_local_stars,
    std::vector<std::vector<int>>& selected_by_chip);
// ==========================================
// Function: Compute one analytic leave-one-out residual and prediction
// Method: Divide the ordinary residual by 1-h after finite leverage guards.
// ==========================================
bool computeAnalyticLOO(
    double observed,
    double fitted,
    double leverage,
    double minimum_denominator,
    double& loo_residual,
    double& loo_model);

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_STAR_SELECTION_HPP
