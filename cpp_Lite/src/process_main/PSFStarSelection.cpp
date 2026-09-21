#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace PSFModel {
namespace Internal {
namespace {
// ==========================================
// Class: Maintain connected components over a compact active-index domain
// Method: Apply path compression and union by rank.
// ==========================================
class DisjointSet {
public:
    // ==========================================
    // Function: Initialize singleton disjoint-set components
    // Method: Assign every compact active index as its own parent with zero rank.
    // ==========================================
    explicit DisjointSet(int size)
        : parent_(static_cast<std::size_t>(size)),
          rank_(static_cast<std::size_t>(size), 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    // ==========================================
    // Function: Find one disjoint-set root
    // Method: Compress every traversed parent link recursively.
    // ==========================================
    int find(int value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    // ==========================================
    // Function: Unite two disjoint-set members
    // Method: Attach the lower-rank root and increment rank on a tie.
    // ==========================================
    void unite(int first, int second) {
        int root_first = find(first);
        int root_second = find(second);
        if (root_first == root_second) return;
        if (rank_[root_first] < rank_[root_second]) {
            std::swap(root_first, root_second);
        }
        parent_[root_second] = root_first;
        if (rank_[root_first] == rank_[root_second]) rank_[root_first]++;
    }

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

// ==========================================
// Structure: Represent one compact same-chip F77 graph edge
// Method: Store the two compact candidate indices joined by the chi cut.
// ==========================================
struct GraphEdge {
    int first = -1;
    int second = -1;
};

// ==========================================
// Structure: Accumulate one connected F77 candidate component
// Method: Preserve compact members in their historical first-seen order.
// ==========================================
struct StarGroup {
    std::vector<int> members;
};

}  // namespace

// ==========================================
// Function: Derive the shared PSF chi-window bounds
// Method: Preserve the legacy inclusive n/4-1 through 3n/4-1 definition.
// ==========================================
PSFChiWindow getPSFChiWindow(int n) {
    if (n <= 0) return {};
    PSFChiWindow window;
    window.first = std::max(0, n / 4 - 1);
    window.last = std::min(n - 1, (n / 4) * 3 - 1);
    return window;
}
// ==========================================
// Function: Count exp(-1)-threshold pixels in one square Fourier-power stamp
// Method: Compare every stamp value with the finite central value times exp(-1).
// ==========================================
int countPSFStarArea(
    const std::vector<float>& power,
    int stamp_side) {
    if (stamp_side <= 0
        || power.size() < static_cast<std::size_t>(stamp_side * stamp_side)) {
        return 0;
    }
    const float center = power[static_cast<std::size_t>(stamp_side / 2)
        * static_cast<std::size_t>(stamp_side)
        + static_cast<std::size_t>(stamp_side / 2)];
    if (!std::isfinite(center)) return 0;
    const float threshold = center * std::exp(-1.0f);
    int star_area = 0;
    for (int index = 0; index < stamp_side * stamp_side; ++index) {
        if (std::isfinite(power[static_cast<std::size_t>(index)])
            && power[static_cast<std::size_t>(index)] >= threshold) {
            ++star_area;
        }
    }
    return star_area;
}

// ==========================================
// Function: Convert an exp(-1) star area to the historical PSF FWHM
// Method: Preserve the legacy area-minus-1e-5 formula and explicit pixel scale.
// ==========================================
double fwhmFromStarArea(
    double star_area,
    int stamp_side,
    double pixel_size) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double area = star_area - 1.0e-5;
    if (!std::isfinite(star_area) || area <= 0.0 || stamp_side <= 0
        || !std::isfinite(pixel_size) || pixel_size <= 0.0) {
        return 0.0;
    }
    const double beta = static_cast<double>(stamp_side) / (2.0 * pi)
        / std::sqrt(area / pi);
    return beta * 2.0 * std::sqrt(2.0 * std::log(2.0)) * pixel_size;
}
// ==========================================
// Function: Compute the exact normalized PSF chi distance
// Method: Apply the legacy sqrt(sum squared difference / mean signed flux)
//         directly to two cached central windows.
// ==========================================
float normalizedChiDistance(
    const std::vector<float>& first,
    const std::vector<float>& second) {
    if (first.empty() || first.size() != second.size()) {
        return std::numeric_limits<float>::infinity();
    }
    double flux = 0.0;
    double squared_difference = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        const double first_value = first[index];
        const double second_value = second[index];
        if (!std::isfinite(first_value) || !std::isfinite(second_value)) {
            return std::numeric_limits<float>::infinity();
        }
        flux += 0.5 * (first_value + second_value);
        const double difference = first_value - second_value;
        squared_difference += difference * difference;
    }
    if (!std::isfinite(flux) || flux <= 0.0
        || !std::isfinite(squared_difference)) {
        return std::numeric_limits<float>::infinity();
    }
    const double distance = std::sqrt(squared_difference / flux);
    return std::isfinite(distance)
        ? static_cast<float>(distance)
        : std::numeric_limits<float>::infinity();
}

// ==========================================
// Function: Prepare one candidate for the compact F77 selection population
// Method: Reject only unsafe inputs, allow a negative finite nonzero full sum,
//         normalize in place, and verify the resulting window remains finite.
// ==========================================
bool prepareF77PSFCandidate(
    int star_index,
    double selection_flag,
    double full_power_sum,
    double size,
    std::vector<float>& chi_window,
    F77PSFCandidateView& candidate) {
    candidate = {};
    if (star_index < 0 || !(selection_flag > 0.0)
        || !std::isfinite(full_power_sum) || full_power_sum == 0.0
        || !std::isfinite(size) || chi_window.empty()
        || !std::all_of(
            chi_window.begin(), chi_window.end(),
            [](float value) { return std::isfinite(value); })) {
        return false;
    }
    const double inverse_sum = 1.0 / full_power_sum;
    for (float& value : chi_window) {
        value = static_cast<float>(static_cast<double>(value) * inverse_sum);
    }
    if (!std::all_of(
            chi_window.begin(), chi_window.end(),
            [](float value) { return std::isfinite(value); })) {
        return false;
    }
    candidate.star_index = star_index;
    candidate.size = size;
    candidate.chi_window = &chi_window;
    return true;
}

// ==========================================
// Function: Test the historical F77 exposure candidate minimum
// Method: Compare only the compact safe population with twice the configured
//         per-exposure minimum.
// ==========================================
bool hasMinimumF77PSFCandidates(
    std::size_t safe_candidate_count,
    int minimum_stars) {
    if (minimum_stars <= 0) return false;
    return safe_candidate_count
        >= 2U * static_cast<std::size_t>(minimum_stars);
}

// ==========================================
// Function: Compute the exact F77 exposure-size and same-chip pair statistics
// Method: Use the one-based floor(2N/3) rank, visit every unordered pair once,
//         and retain threshold pairs only when both sizes meet the rank cut.
// ==========================================
bool computeF77PSFPairStatistics(
    const std::vector<std::vector<F77PSFCandidateView>>& candidates_by_chip,
    F77PSFPairStatistics& statistics) {
    statistics = {};
    std::vector<double> sizes;
    std::size_t total = 0;
    for (const auto& chip : candidates_by_chip) total += chip.size();
    sizes.reserve(total);
    statistics.min_chi.resize(candidates_by_chip.size());
    for (std::size_t chip_index = 0;
         chip_index < candidates_by_chip.size(); ++chip_index) {
        const auto& chip = candidates_by_chip[chip_index];
        statistics.min_chi[chip_index].assign(chip.size(), 1000.0f);
        for (const F77PSFCandidateView& candidate : chip) {
            if (candidate.star_index < 0 || !std::isfinite(candidate.size)
                || candidate.chi_window == nullptr
                || candidate.chi_window->empty()) {
                return false;
            }
            sizes.push_back(candidate.size);
        }
    }
    const std::size_t rank = total * 2U / 3U;
    if (rank == 0U || rank > sizes.size()) return false;
    std::sort(sizes.begin(), sizes.end());
    statistics.size_threshold = sizes[rank - 1U];

    for (std::size_t chip_index = 0;
         chip_index < candidates_by_chip.size(); ++chip_index) {
        const auto& chip = candidates_by_chip[chip_index];
        for (std::size_t first = 0; first + 1U < chip.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < chip.size(); ++second) {
                const float chi = normalizedChiDistance(
                    *chip[first].chi_window, *chip[second].chi_window);
                if (!std::isfinite(chi)) return false;
                statistics.min_chi[chip_index][first] = std::min(
                    statistics.min_chi[chip_index][first], chi);
                statistics.min_chi[chip_index][second] = std::min(
                    statistics.min_chi[chip_index][second], chi);
                if (chip[first].size >= statistics.size_threshold
                    && chip[second].size >= statistics.size_threshold) {
                    statistics.threshold_pair_chi.push_back(chi);
                }
            }
        }
    }
    statistics.valid = true;
    return true;
}
// ==========================================
// Function: Convert one same-chip graph into connected components
// Method: Use disjoint sets over active compact candidate indices.
// ==========================================
std::vector<StarGroup> buildConnectedGroups(
    const std::vector<int>& active_indices,
    const std::vector<GraphEdge>& edges) {
    std::vector<StarGroup> groups;
    if (active_indices.empty()) return groups;

    const int maximum_index = *std::max_element(
        active_indices.begin(), active_indices.end());
    std::vector<int> local_index(static_cast<std::size_t>(maximum_index + 1), -1);
    for (int local = 0; local < static_cast<int>(active_indices.size()); ++local) {
        const int original = active_indices[local];
        if (original >= 0) local_index[original] = local;
    }
    DisjointSet disjoint_set(static_cast<int>(active_indices.size()));
    for (const GraphEdge& edge : edges) {
        if (edge.first < 0 || edge.second < 0
            || edge.first > maximum_index || edge.second > maximum_index) {
            continue;
        }
        const int first_local = local_index[edge.first];
        const int second_local = local_index[edge.second];
        if (first_local >= 0 && second_local >= 0) {
            disjoint_set.unite(first_local, second_local);
        }
    }

    std::vector<int> root_to_group(active_indices.size(), -1);
    for (int local = 0; local < static_cast<int>(active_indices.size()); ++local) {
        const int root = disjoint_set.find(local);
        if (root_to_group[root] < 0) {
            root_to_group[root] = static_cast<int>(groups.size());
            groups.push_back({});
        }
        StarGroup& group = groups[root_to_group[root]];
        const int original = active_indices[local];
        group.members.push_back(original);
    }
    return groups;
}

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
    std::vector<std::vector<int>>& selected_by_chip) {
    selected_by_chip.assign(candidates_by_chip.size(), {});
    if (!statistics.valid || !std::isfinite(chi_threshold)
        || minimum_local_stars <= 0
        || statistics.min_chi.size() != candidates_by_chip.size()) {
        return false;
    }
    for (std::size_t chip_index = 0;
         chip_index < candidates_by_chip.size(); ++chip_index) {
        const auto& chip = candidates_by_chip[chip_index];
        if (statistics.min_chi[chip_index].size() != chip.size()) return false;
        std::vector<int> active;
        for (std::size_t candidate = 0; candidate < chip.size(); ++candidate) {
            if (statistics.min_chi[chip_index][candidate] <= chi_threshold) {
                active.push_back(static_cast<int>(candidate));
            }
        }
        if (static_cast<int>(active.size()) < minimum_local_stars) continue;

        std::vector<GraphEdge> edges;
        for (std::size_t first = 0; first + 1U < active.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < active.size(); ++second) {
                const float chi = normalizedChiDistance(
                    *chip[active[first]].chi_window,
                    *chip[active[second]].chi_window);
                if (!std::isfinite(chi)) return false;
                if (chi <= chi_threshold) {
                    edges.push_back({active[first], active[second]});
                }
            }
        }
        const std::vector<StarGroup> groups = buildConnectedGroups(
            active, edges);
        if (groups.empty()) continue;
        std::size_t largest = 0U;
        for (std::size_t group = 1U; group < groups.size(); ++group) {
            if (groups[group].members.size() > groups[largest].members.size()) {
                largest = group;
            }
        }
        if (static_cast<int>(groups[largest].members.size())
            >= minimum_local_stars) {
            std::vector<int>& selected = selected_by_chip[chip_index];
            selected.reserve(groups[largest].members.size());
            for (int compact_index : groups[largest].members) {
                if (compact_index < 0
                    || compact_index >= static_cast<int>(chip.size())
                    || chip[compact_index].star_index < 0) {
                    return false;
                }
                selected.push_back(chip[compact_index].star_index);
            }
        }
    }
    return true;
}
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
    double& loo_model) {
    constexpr double leverage_tolerance = 1.0e-10;
    if (!std::isfinite(observed) || !std::isfinite(fitted)
        || !std::isfinite(leverage) || !std::isfinite(minimum_denominator)
        || minimum_denominator <= 0.0 || leverage < -leverage_tolerance
        || leverage >= 1.0 - minimum_denominator) {
        return false;
    }
    if (leverage < 0.0) leverage = 0.0;
    const double denominator = 1.0 - leverage;
    loo_residual = (observed - fitted) / denominator;
    loo_model = observed - loo_residual;
    return std::isfinite(loo_residual) && std::isfinite(loo_model);
}

}  // namespace Internal
}  // namespace PSFModel
