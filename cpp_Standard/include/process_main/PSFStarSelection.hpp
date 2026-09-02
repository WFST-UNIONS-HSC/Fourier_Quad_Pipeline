#ifndef PSF_STAR_SELECTION_HPP
#define PSF_STAR_SELECTION_HPP

#include <array>
#include <istream>
#include <string>
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
// Structure: Supply one quality-valid integer star area to locus estimation
// Method: Carry the exact exp(-1) Fourier-pixel count and Gaia support label.
// ==========================================
struct PSFCountSample {
    int star_area = 0;
    bool gaia_matched = false;
};

// ==========================================
// Structure: Configure exposure-wide integer star-area locus estimation
// Method: Group pilot clipping, zero-MAD quantiles, and final-cut controls.
// ==========================================
struct PSFCountLocusConfig {
    double pilot_clip_sigma = 0.0;
    int pilot_clip_iterations = 0;
    double zero_mad_quantile = -1.0;
    double histogram_range_sigma = 0.0;
    double locus_sigma = 0.0;
    int minimum_samples = 0;
    int minimum_gaia_matches = 0;
};

// ==========================================
// Structure: Publish a robust exposure-wide integer star-area locus
// Method: Store the median center, side-specific widths, and strict-cut bounds.
// ==========================================
struct PSFCountLocus {
    bool valid = false;
    double center = 0.0;
    double lower_width = 0.0;
    double upper_width = 0.0;
    double lower = 0.0;
    double upper = 0.0;
};

// ==========================================
// Structure: Publish exact integer star-area locus diagnostics
// Method: Keep raw counts immutable, expose the hole-filled working histogram,
//         and retain pilot, range, Gaia, and selected-peak accounting.
// ==========================================
struct PSFCountLocusDiagnostics {
    int sample_count = 0;
    int gaia_match_count = 0;
    bool pilot_uses_gaia = false;
    int pilot_input_count = 0;
    int pilot_retained_count = 0;
    double pilot_center = 0.0;
    double pilot_width = 0.0;
    double pilot_lower = 0.0;
    double pilot_upper = 0.0;
    bool pilot_uses_quantile_range = false;
    bool pilot_rejected_zero_mad_clip = false;
    int histogram_sample_count = 0;
    int histogram_below_count = 0;
    int histogram_above_count = 0;
    int gaia_histogram_sample_count = 0;
    int gaia_histogram_below_count = 0;
    int gaia_histogram_above_count = 0;
    int selected_group_count = 0;
    bool has_gaia_median = false;
    double gaia_median = 0.0;
    int histogram_first_count = 0;
    int peak_bin = -1;
    std::vector<double> histogram;
    std::vector<double> working_histogram;
    std::vector<double> smoothed_histogram;
    std::vector<double> gaia_histogram;
    std::vector<double> selected_group_histogram;
};

// ==========================================
// Function: Select one integer star-area histogram peak deterministically
// Method: Anchor Gaia eligibility to the global nearest distance plus one count,
//         then rank by Gaia count, density, exact distance, and lower count.
// ==========================================
int selectPSFCountPeak(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    const std::vector<double>& gaia_histogram,
    int histogram_first_count,
    double pilot_center,
    bool pilot_uses_gaia);

// ==========================================
// Function: Fill bounded one- or two-level holes in a working histogram
// Method: Detect zero runs only in the immutable raw histogram and linearly
//         interpolate between positive endpoints without chaining.
// ==========================================
std::vector<double> interpolateShortInternalHoles(
    const std::vector<double>& histogram);

// ==========================================
// Function: Estimate an exposure-wide integer star-area locus
// Method: Build one-count bins, fill only short internal holes for smoothing,
//         select a count peak, and refine its basin with asymmetric count MAD.
// ==========================================
bool estimatePSFCountLocus(
    const std::vector<PSFCountSample>& samples,
    const PSFCountLocusConfig& config,
    PSFCountLocus& locus,
    PSFCountLocusDiagnostics* diagnostics = nullptr);

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

// ==========================================
// Function: Populate the shared-group integer-area histogram
// Method: Count all selected stars and bin only values on the science grid
//         without changing upstream count-locus diagnostics.
// ==========================================
void populateSelectedGroupCountHistogram(
    const std::vector<int>& selected_star_areas,
    PSFCountLocusDiagnostics& diagnostics);

enum class AstrometryGaiaReadStatus {
    Accepted,
    Empty,
    Malformed
};

// ==========================================
// Function: Parse matched Gaia image positions from one astro stream
// Method: Read the two WCS lines, the matched/user/reference counts, and exactly
//         n_matched finite RA/Dec/x/y rows while retaining only x/y.
// ==========================================
AstrometryGaiaReadStatus parseAstrometryGaiaPositions(
    std::istream& input,
    std::vector<std::array<double, 2>>& gaia_xy,
    std::string& error);

// ==========================================
// Function: Match one candidate to the nearest same-chip Gaia position
// Method: Accept when the minimum squared image-plane distance is within the
//         configured radius; do not impose one-to-one assignment.
// ==========================================
bool hasNearestGaiaMatch(
    double x,
    double y,
    const std::vector<std::array<double, 2>>& gaia_xy,
    double radius_pixels);

// ==========================================
// Structure: Store one exact same-chip nearest-neighbour relation
// Method: Retain the original candidate index and its normalized Fourier chi.
// ==========================================
struct NeighborEdge {
    int star_index = -1;
    float chi = 0.0f;
};

// ==========================================
// Structure: Identify one exposure-level minChi reference candidate
// Method: Carry deterministic chip/star keys and the legacy index-7 size rank.
// ==========================================
struct MinChiReferenceCandidate {
    int chip_index = -1;
    int star_index = -1;
    double size = 0.0;
};

// ==========================================
// Function: Select capped exposure-wide large-size minChi references
// Method: Sort the top fraction by size/chip/star, then apply a per-chip cap
//         without filling from candidates below the exposure-wide pool.
// ==========================================
std::vector<MinChiReferenceCandidate> selectMinChiReferenceStars(
    const std::vector<MinChiReferenceCandidate>& locus_candidates,
    double reference_fraction,
    int maximum_per_chip);

// ==========================================
// Structure: View one candidate during the same-chip minChi pair pass
// Method: Borrow its normalized chi window and carry locus/reference labels.
// ==========================================
struct MinChiCandidateView {
    const std::vector<float>* chi_window = nullptr;
    bool in_size_locus = false;
    bool is_reference = false;
};

// ==========================================
// Structure: Return nearest-neighbour distances and reference-pair samples
// Method: Align minChi with input candidates and retain every qualifying
//         unordered pair distance exactly once for exposure thresholding.
// ==========================================
struct MinChiPairResult {
    std::vector<float> min_chi;
    std::vector<float> threshold_pair_chi;
};

// ==========================================
// Function: Compute one chip's minChi values and threshold-pair sample
// Method: Visit every unordered locus-locus pair once, update both endpoints,
//         and sample the distance when either endpoint is a reference.
// ==========================================
MinChiPairResult computeMinChiAndThresholdPairs(
    const std::vector<MinChiCandidateView>& candidates);

// ==========================================
// Function: Compute the exact normalized PSF chi distance
// Method: Apply the legacy sqrt(sum squared difference / mean signed flux)
//         directly to two cached central windows.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second);

// ==========================================
// Function: Maintain one exact sorted top-K neighbour list
// Method: Insert or improve the candidate edge, sort by chi/index, and truncate.
// ==========================================
void updateTopK(
    std::vector<NeighborEdge>& neighbours,
    int star_index,
    float chi,
    int k);

// ==========================================
// Function: Rebuild exact top-K neighbours on one active survivor set
// Method: Clear every stale list, compare only active same-chip cached windows,
//         and refill bounded neighbours in exact distance/index order.
// ==========================================
template <typename CandidateSelectionState>
void rebuildActiveKNN(
    const std::vector<int>& active_indices,
    std::vector<CandidateSelectionState>& selection,
    int k) {
    for (CandidateSelectionState& candidate : selection) {
        candidate.knn.clear();
    }
    if (k <= 0) return;

    for (std::size_t first = 0; first + 1 < active_indices.size(); ++first) {
        const int first_index = active_indices[first];
        if (first_index < 0
            || first_index >= static_cast<int>(selection.size())) {
            continue;
        }
        for (std::size_t second = first + 1;
             second < active_indices.size(); ++second) {
            const int second_index = active_indices[second];
            if (second_index < 0 || second_index == first_index
                || second_index >= static_cast<int>(selection.size())) {
                continue;
            }
            const float chi = normalizedChiDistance(
                selection[first_index].chi_window,
                selection[second_index].chi_window);
            updateTopK(
                selection[first_index].knn, second_index, chi, k);
            updateTopK(
                selection[second_index].knn, first_index, chi, k);
        }
    }
}

// ==========================================
// Structure: Describe one undirected grouping-graph edge
// Method: Address both endpoints by their original per-chip candidate indices.
// ==========================================
struct GraphEdge {
    int first = -1;
    int second = -1;
};

// ==========================================
// Structure: Return one connected stellar component
// Method: Preserve original candidate indices and count Gaia supporting labels.
// ==========================================
struct StarGroup {
    std::vector<int> members;
    int gaia_count = 0;
};

// ==========================================
// Function: Extract mutual-KNN graph edges among active candidates
// Method: Keep an undirected edge only when both retained top-K lists contain
//         the opposite endpoint and both endpoints survived the shared cut.
// ==========================================
std::vector<GraphEdge> buildMutualKNNEdges(
    const std::vector<int>& active_indices,
    const std::vector<std::vector<NeighborEdge>>& neighbours_by_star);

// ==========================================
// Function: Convert one same-chip graph into connected components
// Method: Use disjoint sets over active original indices and attach Gaia counts.
// ==========================================
std::vector<StarGroup> buildConnectedGroups(
    const std::vector<int>& active_indices,
    const std::vector<GraphEdge>& edges,
    const std::vector<bool>& gaia_matched);

// ==========================================
// Function: Select the shared main and eligible secondary stellar groups
// Method: Always keep the largest component and keep another component only
//         when both its relative size and Gaia-count requirements pass.
// ==========================================
std::vector<int> selectMainAndSecondaryGroups(
    const std::vector<StarGroup>& groups,
    double minimum_size_ratio,
    int minimum_gaia_count);

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

// ==========================================
// Function: Convert raw analytic PRESS to its leverage-standardized score
// Method: Multiply by sqrt(1-h) after the same finite denominator guard used
//         by the analytic leave-one-out calculation.
// ==========================================
bool computeLeverageStandardizedPress(
    double raw_press,
    double leverage,
    double minimum_denominator,
    double& standardized_press);

enum class PressRemovalDecision {
    Disabled,
    NoOutliers,
    TooManyOutliers,
    WouldUnderrunMinimum,
    Apply
};

// ==========================================
// Function: Decide whether an optional PRESS removal may be attempted
// Method: Apply the rejection switch, outlier count, removal cap, and minimum
//         retained-star guard without mutating any fitting or selection state.
// ==========================================
PressRemovalDecision decidePressRemoval(
    bool rejection_enabled,
    int initial_count,
    int flagged_count,
    int minimum_count,
    int maximum_removals);

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_STAR_SELECTION_HPP
