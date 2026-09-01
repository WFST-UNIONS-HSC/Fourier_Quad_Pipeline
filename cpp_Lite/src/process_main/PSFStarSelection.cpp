#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace PSFModel {
namespace Internal {
namespace {

// ==========================================
// Function: Return the median of a sorted finite vector
// Method: Select the middle element or average the two central elements.
// ==========================================
double sortedMedian(const std::vector<double>& sorted) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t middle = sorted.size() / 2;
    if (sorted.size() % 2 == 1) return sorted[middle];
    return 0.5 * (sorted[middle - 1] + sorted[middle]);
}

// ==========================================
// Function: Interpolate one quantile from a sorted finite vector
// Method: Map the bounded fraction onto [0, N-1] and linearly blend neighbours.
// ==========================================
double sortedQuantile(
    const std::vector<double>& sorted,
    double fraction) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    fraction = std::clamp(fraction, 0.0, 1.0);
    const double position = fraction
        * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    if (lower == upper) return sorted[lower];
    const double weight = position - static_cast<double>(lower);
    return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
}

// ==========================================
// Function: Estimate median and scaled median absolute deviation
// Method: Sort finite samples and their deviations and multiply MAD by 1.4826.
// ==========================================
bool medianAndMad(
    const std::vector<double>& values,
    double& median,
    double& scaled_mad) {
    if (values.empty()) return false;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    median = sortedMedian(sorted);
    std::vector<double> deviations;
    deviations.reserve(sorted.size());
    for (double value : sorted) deviations.push_back(std::abs(value - median));
    std::sort(deviations.begin(), deviations.end());
    scaled_mad = 1.4826 * sortedMedian(deviations);
    return std::isfinite(median) && std::isfinite(scaled_mad);
}

// ==========================================
// Function: Estimate median and side-specific scaled median deviations
// Method: Exclude center duplicates from both side MADs and apply the fixed
//         one-count measurement-resolution floor independently to both sides.
// ==========================================
bool medianAndAsymmetricMad(
    const std::vector<double>& values,
    double& median,
    double& lower_width,
    double& upper_width) {
    if (values.empty()) return false;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    median = sortedMedian(sorted);

    std::vector<double> lower_deviations;
    std::vector<double> upper_deviations;
    lower_deviations.reserve(sorted.size());
    upper_deviations.reserve(sorted.size());
    for (double value : sorted) {
        if (value < median) {
            lower_deviations.push_back(median - value);
        } else if (value > median) {
            upper_deviations.push_back(value - median);
        }
    }
    std::sort(lower_deviations.begin(), lower_deviations.end());
    std::sort(upper_deviations.begin(), upper_deviations.end());
    lower_width = lower_deviations.empty()
        ? 0.0
        : 1.4826 * sortedMedian(lower_deviations);
    upper_width = upper_deviations.empty()
        ? 0.0
        : 1.4826 * sortedMedian(upper_deviations);
    lower_width = std::max(lower_width, 1.0);
    upper_width = std::max(upper_width, 1.0);
    return std::isfinite(median)
        && std::isfinite(lower_width) && lower_width > 0.0
        && std::isfinite(upper_width) && upper_width > 0.0;
}

// ==========================================
// Structure: Publish one first-layer integer star-area pilot and its bounds
// Method: Keep center/width, retained count, and zero-MAD branch decisions atomic.
// ==========================================
struct PSFCountPilotEstimate {
    double center = 0.0;
    double width = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    int retained_count = 0;
    bool uses_quantile_range = false;
    bool rejected_zero_mad_clip = false;
};

// ==========================================
// Function: Estimate a quantization-safe first-layer integer star-area pilot
// Method: Use configurable symmetric Q(q)--Q(1-q) bounds for an initially zero
//         MAD; otherwise accept clips only while their next MAD remains positive.
// ==========================================
bool estimatePSFCountPilot(
    const std::vector<double>& sorted_input,
    double clip_sigma,
    int max_iterations,
    double zero_mad_quantile,
    double range_sigma,
    PSFCountPilotEstimate& pilot) {
    pilot = {};
    if (sorted_input.empty() || !std::isfinite(clip_sigma)
        || clip_sigma <= 0.0 || max_iterations <= 0
        || !std::isfinite(zero_mad_quantile)
        || zero_mad_quantile < 0.0 || zero_mad_quantile >= 0.5
        || !std::isfinite(range_sigma) || range_sigma <= 0.0) {
        return false;
    }

    std::vector<double> population = sorted_input;
    double center = 0.0;
    double width = 0.0;
    if (!medianAndMad(population, center, width)) return false;

    if (width <= 0.0) {
        pilot.center = center;
        pilot.width = 0.0;
        pilot.lower = sortedQuantile(population, zero_mad_quantile);
        pilot.upper = sortedQuantile(population, 1.0 - zero_mad_quantile);
        pilot.retained_count = static_cast<int>(population.size());
        pilot.uses_quantile_range = true;
        return std::isfinite(pilot.center)
            && std::isfinite(pilot.width)
            && pilot.width >= 0.0
            && std::isfinite(pilot.lower)
            && std::isfinite(pilot.upper)
            && pilot.upper >= pilot.lower
            && pilot.retained_count > 0;
    }

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (!medianAndMad(population, center, width)) return false;
        std::vector<double> clipped;
        clipped.reserve(population.size());
        for (double value : population) {
            if (std::abs(value - center) <= clip_sigma * width) {
                clipped.push_back(value);
            }
        }
        if (clipped.empty()) return false;
        if (clipped.size() == population.size()) break;

        double next_center = 0.0;
        double next_width = 0.0;
        if (!medianAndMad(clipped, next_center, next_width)) return false;
        if (next_width <= 0.0) {
            pilot.rejected_zero_mad_clip = true;
            break;
        }
        population.swap(clipped);
    }

    if (!medianAndMad(population, pilot.center, pilot.width)) return false;
    pilot.retained_count = static_cast<int>(population.size());
    if (pilot.width <= 0.0) {
        pilot.uses_quantile_range = true;
        pilot.lower = sortedQuantile(population, zero_mad_quantile);
        pilot.upper = sortedQuantile(population, 1.0 - zero_mad_quantile);
    } else {
        pilot.lower = pilot.center - range_sigma * pilot.width;
        pilot.upper = pilot.center + range_sigma * pilot.width;
    }
    return std::isfinite(pilot.center)
        && std::isfinite(pilot.width)
        && pilot.width >= 0.0
        && std::isfinite(pilot.lower)
        && std::isfinite(pilot.upper)
        && pilot.upper >= pilot.lower
        && pilot.retained_count > 0;
}

// ==========================================
// Function: Test whether a top-K list contains one candidate index
// Method: Perform a bounded linear scan over the short neighbour vector.
// ==========================================
bool containsNeighbour(
    const std::vector<NeighborEdge>& neighbours,
    int star_index) {
    return std::any_of(
        neighbours.begin(), neighbours.end(),
        [star_index](const NeighborEdge& edge) {
            return edge.star_index == star_index;
        });
}

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
// Function: Select one integer star-area histogram peak deterministically
// Method: Use density without Gaia; with Gaia, anchor eligibility to the global
//         nearest distance plus one count and apply the deterministic ranking.
// ==========================================
int selectPSFCountPeak(
    const std::vector<int>& peaks,
    const std::vector<double>& smoothed_histogram,
    const std::vector<double>& gaia_histogram,
    int histogram_first_count,
    double pilot_center,
    bool pilot_uses_gaia) {
    if (peaks.empty() || smoothed_histogram.empty()
        || gaia_histogram.size() != smoothed_histogram.size()
        || histogram_first_count <= 0
        || !std::isfinite(pilot_center)) {
        return -1;
    }
    for (int peak : peaks) {
        if (peak < 0
            || peak >= static_cast<int>(smoothed_histogram.size())) {
            return -1;
        }
    }

    int selected = peaks.front();
    if (!pilot_uses_gaia) {
        for (int candidate : peaks) {
            if (smoothed_histogram[candidate]
                > smoothed_histogram[selected]) {
                selected = candidate;
            }
        }
        return selected;
    }

    double minimum_distance = std::numeric_limits<double>::infinity();
    for (int peak : peaks) {
        const double distance = std::abs(
            static_cast<double>(histogram_first_count + peak) - pilot_center);
        minimum_distance = std::min(minimum_distance, distance);
    }

    selected = -1;
    double selected_distance = std::numeric_limits<double>::infinity();
    for (int candidate : peaks) {
        const double candidate_distance = std::abs(
            static_cast<double>(histogram_first_count + candidate)
                - pilot_center);
        if (candidate_distance > minimum_distance + 1.0) {
            continue;
        }
        if (selected < 0) {
            selected = candidate;
            selected_distance = candidate_distance;
            continue;
        }

        bool replace = false;
        if (gaia_histogram[candidate] != gaia_histogram[selected]) {
            replace = gaia_histogram[candidate]
                > gaia_histogram[selected];
        } else if (smoothed_histogram[candidate]
                   != smoothed_histogram[selected]) {
            replace = smoothed_histogram[candidate]
                > smoothed_histogram[selected];
        } else if (candidate_distance != selected_distance) {
            replace = candidate_distance < selected_distance;
        } else {
            replace = candidate < selected;
        }
        if (replace) {
            selected = candidate;
            selected_distance = candidate_distance;
        }
    }
    return selected;
}

// ==========================================
// Function: Fill bounded one- or two-level holes in a working histogram
// Method: Detect every zero run from the immutable raw input, then interpolate
//         only short internal runs between positive raw endpoints.
// ==========================================
std::vector<double> interpolateShortInternalHoles(
    const std::vector<double>& histogram) {
    std::vector<double> working = histogram;
    std::size_t index = 0;
    while (index < histogram.size()) {
        if (histogram[index] != 0.0) {
            ++index;
            continue;
        }
        const std::size_t first = index;
        while (index < histogram.size() && histogram[index] == 0.0) ++index;
        const std::size_t length = index - first;
        if (first == 0 || index >= histogram.size() || length > 2
            || histogram[first - 1] <= 0.0 || histogram[index] <= 0.0) {
            continue;
        }
        const double left = histogram[first - 1];
        const double right = histogram[index];
        for (std::size_t offset = 0; offset < length; ++offset) {
            const double fraction = static_cast<double>(offset + 1)
                / static_cast<double>(length + 1);
            working[first + offset] = left + fraction * (right - left);
        }
    }
    return working;
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
// Function: Estimate an exposure-wide integer star-area locus
// Method: Build one-count bins, fill only short internal holes for smoothing,
//         select a count peak, and refine its basin with asymmetric count MAD.
// ==========================================
bool estimatePSFCountLocus(
    const std::vector<PSFCountSample>& samples,
    const PSFCountLocusConfig& config,
    PSFCountLocus& locus,
    PSFCountLocusDiagnostics* diagnostics) {
    locus = {};
    if (diagnostics != nullptr) *diagnostics = {};
    if (!std::isfinite(config.pilot_clip_sigma)
        || config.pilot_clip_sigma <= 0.0
        || config.pilot_clip_iterations <= 0
        || !std::isfinite(config.zero_mad_quantile)
        || config.zero_mad_quantile < 0.0
        || config.zero_mad_quantile >= 0.5
        || !std::isfinite(config.histogram_range_sigma)
        || config.histogram_range_sigma <= 0.0
        || !std::isfinite(config.locus_sigma)
        || config.locus_sigma <= 0.0
        || config.minimum_samples <= 0
        || config.minimum_gaia_matches <= 0) {
        return false;
    }

    std::vector<double> values;
    std::vector<double> gaia_values;
    values.reserve(samples.size());
    gaia_values.reserve(samples.size());
    for (const PSFCountSample& sample : samples) {
        if (sample.star_area <= 0) continue;
        values.push_back(static_cast<double>(sample.star_area));
        if (sample.gaia_matched) {
            gaia_values.push_back(static_cast<double>(sample.star_area));
        }
    }
    if (diagnostics != nullptr) {
        diagnostics->sample_count = static_cast<int>(values.size());
        diagnostics->gaia_match_count = static_cast<int>(gaia_values.size());
    }
    if (static_cast<int>(values.size()) < config.minimum_samples) return false;

    std::sort(values.begin(), values.end());
    std::sort(gaia_values.begin(), gaia_values.end());
    const bool has_gaia_support =
        static_cast<int>(gaia_values.size()) >= config.minimum_gaia_matches;
    if (diagnostics != nullptr && has_gaia_support) {
        diagnostics->has_gaia_median = true;
        diagnostics->gaia_median = sortedMedian(gaia_values);
    }

    bool pilot_uses_gaia = has_gaia_support;
    int pilot_input_count = pilot_uses_gaia
        ? static_cast<int>(gaia_values.size())
        : static_cast<int>(values.size());
    PSFCountPilotEstimate pilot;
    bool pilot_valid = estimatePSFCountPilot(
        pilot_uses_gaia ? gaia_values : values,
        config.pilot_clip_sigma,
        config.pilot_clip_iterations,
        config.zero_mad_quantile,
        config.histogram_range_sigma,
        pilot);
    if (pilot_uses_gaia
        && (!pilot_valid
            || pilot.retained_count < config.minimum_gaia_matches)) {
        pilot_uses_gaia = false;
        pilot_input_count = static_cast<int>(values.size());
        pilot = {};
        pilot_valid = estimatePSFCountPilot(
            values,
            config.pilot_clip_sigma,
            config.pilot_clip_iterations,
            config.zero_mad_quantile,
            config.histogram_range_sigma,
            pilot);
    }
    if (!pilot_valid) return false;

    if (diagnostics != nullptr) {
        diagnostics->pilot_uses_gaia = pilot_uses_gaia;
        diagnostics->pilot_input_count = pilot_input_count;
        diagnostics->pilot_retained_count = pilot.retained_count;
        diagnostics->pilot_center = pilot.center;
        diagnostics->pilot_width = pilot.width;
        diagnostics->pilot_lower = pilot.lower;
        diagnostics->pilot_upper = pilot.upper;
        diagnostics->pilot_uses_quantile_range = pilot.uses_quantile_range;
        diagnostics->pilot_rejected_zero_mad_clip =
            pilot.rejected_zero_mad_clip;
    }

    int histogram_first_count = std::max(
        1, static_cast<int>(std::ceil(pilot.lower)));
    int histogram_last_count = static_cast<int>(std::floor(pilot.upper));
    if (histogram_first_count > histogram_last_count) {
        const std::vector<double>& source = pilot_uses_gaia
            ? gaia_values : values;
        if (source.empty()) return false;
        int nearest = static_cast<int>(std::llround(source.front()));
        double nearest_distance = std::abs(
            static_cast<double>(nearest) - pilot.center);
        for (double value : source) {
            const int observed = static_cast<int>(std::llround(value));
            const double distance = std::abs(value - pilot.center);
            if (distance < nearest_distance
                || (distance == nearest_distance && observed < nearest)) {
                nearest = observed;
                nearest_distance = distance;
            }
        }
        histogram_first_count = nearest;
        histogram_last_count = nearest;
    }
    const int histogram_bin_count =
        histogram_last_count - histogram_first_count + 1;
    if (histogram_bin_count <= 0) return false;

    std::vector<double> histogram(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    int histogram_sample_count = 0;
    int histogram_below_count = 0;
    int histogram_above_count = 0;
    for (double value : values) {
        const int star_area = static_cast<int>(std::llround(value));
        if (star_area < histogram_first_count) {
            histogram_below_count++;
            continue;
        }
        if (star_area > histogram_last_count) {
            histogram_above_count++;
            continue;
        }
        const int bin = star_area - histogram_first_count;
        histogram[bin] += 1.0;
        histogram_sample_count++;
    }

    std::vector<double> gaia_histogram(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    int gaia_histogram_sample_count = 0;
    int gaia_histogram_below_count = 0;
    int gaia_histogram_above_count = 0;
    for (double value : gaia_values) {
        const int star_area = static_cast<int>(std::llround(value));
        if (star_area < histogram_first_count) {
            gaia_histogram_below_count++;
            continue;
        }
        if (star_area > histogram_last_count) {
            gaia_histogram_above_count++;
            continue;
        }
        const int bin = star_area - histogram_first_count;
        gaia_histogram[bin] += 1.0;
        gaia_histogram_sample_count++;
    }

    const std::vector<double> working =
        interpolateShortInternalHoles(histogram);
    const int offsets[] = {-2, -1, 0, 1, 2};
    const double weights[] = {1.0, 2.0, 3.0, 2.0, 1.0};
    std::vector<double> smoothed(
        static_cast<std::size_t>(histogram_bin_count), 0.0);
    for (int bin = 0; bin < histogram_bin_count; ++bin) {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (int offset_index = 0; offset_index < 5; ++offset_index) {
            const int neighbour = bin + offsets[offset_index];
            if (neighbour < 0 || neighbour >= histogram_bin_count) continue;
            weighted_sum += weights[offset_index] * working[neighbour];
            weight_sum += weights[offset_index];
        }
        smoothed[bin] = weighted_sum / weight_sum;
    }

    std::vector<int> peaks;
    for (int bin = 0; bin < histogram_bin_count; ++bin) {
        const double left = bin == 0 ? -1.0 : smoothed[bin - 1];
        const double right = bin + 1 == histogram_bin_count
            ? -1.0
            : smoothed[bin + 1];
        if (smoothed[bin] > 0.0 && smoothed[bin] >= left && smoothed[bin] >= right) {
            peaks.push_back(bin);
        }
    }
    if (peaks.empty()) return false;

    const int peak_bin = selectPSFCountPeak(
        peaks,
        smoothed,
        gaia_histogram,
        histogram_first_count,
        pilot.center,
        pilot_uses_gaia);
    if (peak_bin < 0) return false;
    if (diagnostics != nullptr) {
        diagnostics->histogram_sample_count = histogram_sample_count;
        diagnostics->histogram_below_count = histogram_below_count;
        diagnostics->histogram_above_count = histogram_above_count;
        diagnostics->gaia_histogram_sample_count =
            gaia_histogram_sample_count;
        diagnostics->gaia_histogram_below_count =
            gaia_histogram_below_count;
        diagnostics->gaia_histogram_above_count =
            gaia_histogram_above_count;
        diagnostics->histogram_first_count = histogram_first_count;
        diagnostics->peak_bin = peak_bin;
        diagnostics->histogram = histogram;
        diagnostics->working_histogram = working;
        diagnostics->smoothed_histogram = smoothed;
        diagnostics->gaia_histogram = gaia_histogram;
    }

    int basin_first = peak_bin;
    while (basin_first > 0
           && smoothed[basin_first - 1] <= smoothed[basin_first]) {
        basin_first--;
    }
    int basin_last = peak_bin;
    while (basin_last + 1 < histogram_bin_count
           && smoothed[basin_last + 1] <= smoothed[basin_last]) {
        basin_last++;
    }

    const int basin_low = histogram_first_count + basin_first;
    const int basin_high = histogram_first_count + basin_last;
    std::vector<double> population;
    for (double value : values) {
        if (value >= static_cast<double>(basin_low)
            && value <= static_cast<double>(basin_high)) {
            population.push_back(value);
        }
    }
    if (population.size() < 3) {
        const int fallback_first = std::max(0, peak_bin - 2);
        const int fallback_last = std::min(
            histogram_bin_count - 1, peak_bin + 2);
        const int fallback_low = histogram_first_count + fallback_first;
        const int fallback_high = histogram_first_count + fallback_last;
        population.clear();
        for (double value : values) {
            if (value >= static_cast<double>(fallback_low)
                && value <= static_cast<double>(fallback_high)) {
                population.push_back(value);
            }
        }
    }
    if (population.empty()) return false;

    double center = 0.0;
    double lower_width = 0.0;
    double upper_width = 0.0;
    for (int iteration = 0; iteration < 2; ++iteration) {
        if (!medianAndAsymmetricMad(
                population, center, lower_width, upper_width)) {
            return false;
        }
        const double clip_lower =
            center - config.locus_sigma * lower_width;
        const double clip_upper =
            center + config.locus_sigma * upper_width;
        std::vector<double> clipped;
        clipped.reserve(population.size());
        for (double value : population) {
            if (value >= clip_lower && value <= clip_upper) {
                clipped.push_back(value);
            }
        }
        if (clipped.empty() || clipped.size() == population.size()) break;
        population.swap(clipped);
    }
    if (!medianAndAsymmetricMad(
            population, center, lower_width, upper_width)) {
        return false;
    }

    locus.valid = true;
    locus.center = center;
    locus.lower_width = lower_width;
    locus.upper_width = upper_width;
    locus.lower = center - config.locus_sigma * lower_width;
    locus.upper = center + config.locus_sigma * upper_width;
    return std::isfinite(locus.lower) && std::isfinite(locus.upper);
}

// ==========================================
// Function: Build the unchanged FWHM SVG view from isolated display samples
// Method: Map count science markers to FWHM and independently histogram FWHM.
// ==========================================
bool buildFWHMLocusDisplay(
    const std::vector<FWHMDisplaySample>& samples,
    const PSFCountLocus& count_locus,
    const PSFCountLocusDiagnostics& count_diagnostics,
    double locus_sigma,
    int stamp_side,
    double pixel_size,
    FWHMDisplayLocus& display_locus,
    FWHMDisplayDiagnostics& display_diagnostics) {
    display_locus = {};
    display_diagnostics = {};
    const int bin_count = static_cast<int>(
        count_diagnostics.histogram.size());
    if (!count_locus.valid || bin_count <= 0 || locus_sigma <= 0.0
        || count_diagnostics.peak_bin < 0
        || count_diagnostics.peak_bin >= bin_count) {
        return false;
    }
    const int first_count = count_diagnostics.histogram_first_count;
    const int last_count = first_count + bin_count - 1;
    const auto map_count = [&](double value) {
        return fwhmFromStarArea(
            std::max(value, 0.5), stamp_side, pixel_size);
    };

    display_locus.center = map_count(count_locus.center);
    display_locus.lower = map_count(count_locus.upper);
    display_locus.upper = map_count(count_locus.lower);
    display_locus.lower_width =
        (display_locus.center - display_locus.lower) / locus_sigma;
    display_locus.upper_width =
        (display_locus.upper - display_locus.center) / locus_sigma;
    display_locus.valid = std::isfinite(display_locus.center)
        && display_locus.center > 0.0
        && display_locus.lower < display_locus.upper;
    if (!display_locus.valid) return false;

    display_diagnostics.sample_count = count_diagnostics.sample_count;
    display_diagnostics.gaia_match_count = count_diagnostics.gaia_match_count;
    display_diagnostics.pilot_uses_gaia =
        count_diagnostics.pilot_uses_gaia;
    display_diagnostics.pilot_input_count =
        count_diagnostics.pilot_input_count;
    display_diagnostics.pilot_retained_count =
        count_diagnostics.pilot_retained_count;
    display_diagnostics.pilot_center =
        map_count(count_diagnostics.pilot_center);
    display_diagnostics.pilot_lower =
        map_count(count_diagnostics.pilot_upper);
    display_diagnostics.pilot_upper =
        map_count(count_diagnostics.pilot_lower);
    display_diagnostics.pilot_width = count_diagnostics.pilot_width > 0.0
        ? std::max(
            display_diagnostics.pilot_center
                - map_count(count_diagnostics.pilot_center
                    + count_diagnostics.pilot_width),
            map_count(count_diagnostics.pilot_center
                    - count_diagnostics.pilot_width)
                - display_diagnostics.pilot_center)
        : 0.0;
    display_diagnostics.pilot_uses_quantile_range =
        count_diagnostics.pilot_uses_quantile_range;
    display_diagnostics.pilot_rejected_zero_mad_clip =
        count_diagnostics.pilot_rejected_zero_mad_clip;
    display_diagnostics.has_gaia_median =
        count_diagnostics.has_gaia_median;
    if (display_diagnostics.has_gaia_median) {
        display_diagnostics.gaia_median =
            map_count(count_diagnostics.gaia_median);
    }
    display_diagnostics.histogram_lower =
        map_count(static_cast<double>(last_count) + 0.5);
    display_diagnostics.histogram_upper =
        map_count(static_cast<double>(first_count) - 0.5);
    display_diagnostics.peak_value = map_count(
        static_cast<double>(first_count + count_diagnostics.peak_bin));
    if (!(display_diagnostics.histogram_upper
            > display_diagnostics.histogram_lower)) {
        return false;
    }

    display_diagnostics.histogram.assign(
        static_cast<std::size_t>(bin_count), 0.0);
    display_diagnostics.gaia_histogram.assign(
        static_cast<std::size_t>(bin_count), 0.0);
    const double span = display_diagnostics.histogram_upper
        - display_diagnostics.histogram_lower;
    for (const FWHMDisplaySample& sample : samples) {
        if (!std::isfinite(sample.fwhm) || sample.fwhm <= 0.0) continue;
        int* below = sample.gaia_matched
            ? &display_diagnostics.gaia_histogram_below_count
            : nullptr;
        int* above = sample.gaia_matched
            ? &display_diagnostics.gaia_histogram_above_count
            : nullptr;
        if (sample.fwhm < display_diagnostics.histogram_lower) {
            ++display_diagnostics.histogram_below_count;
            if (below != nullptr) ++(*below);
            continue;
        }
        if (sample.fwhm > display_diagnostics.histogram_upper) {
            ++display_diagnostics.histogram_above_count;
            if (above != nullptr) ++(*above);
            continue;
        }
        int bin = static_cast<int>((sample.fwhm
            - display_diagnostics.histogram_lower) / span * bin_count);
        bin = std::clamp(bin, 0, bin_count - 1);
        display_diagnostics.histogram[static_cast<std::size_t>(bin)] += 1.0;
        ++display_diagnostics.histogram_sample_count;
        if (sample.gaia_matched) {
            display_diagnostics.gaia_histogram[static_cast<std::size_t>(bin)]
                += 1.0;
            ++display_diagnostics.gaia_histogram_sample_count;
        }
    }

    const int offsets[] = {-2, -1, 0, 1, 2};
    const double weights[] = {1.0, 2.0, 3.0, 2.0, 1.0};
    display_diagnostics.smoothed_histogram.assign(
        static_cast<std::size_t>(bin_count), 0.0);
    for (int bin = 0; bin < bin_count; ++bin) {
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (int offset_index = 0; offset_index < 5; ++offset_index) {
            const int neighbour = bin + offsets[offset_index];
            if (neighbour < 0 || neighbour >= bin_count) continue;
            weighted_sum += weights[offset_index]
                * display_diagnostics.histogram[neighbour];
            weight_sum += weights[offset_index];
        }
        display_diagnostics.smoothed_histogram[bin] =
            weighted_sum / weight_sum;
    }
    return true;
}

// ==========================================
// Function: Populate the shared-group FWHM overlay on the display bin grid
// Method: Count every selected value and bin finite in-range values without
//         changing science or upstream display diagnostics.
// ==========================================
void populateSelectedGroupFWHMHistogram(
    const std::vector<double>& selected_fwhm,
    FWHMDisplayDiagnostics& diagnostics) {
    diagnostics.selected_group_count =
        static_cast<int>(selected_fwhm.size());
    diagnostics.selected_group_histogram.assign(
        diagnostics.histogram.size(), 0.0);
    if (diagnostics.histogram.empty()
        || !(diagnostics.histogram_upper > diagnostics.histogram_lower)) {
        return;
    }
    const double bin_width =
        (diagnostics.histogram_upper - diagnostics.histogram_lower)
        / static_cast<double>(diagnostics.histogram.size());
    if (!std::isfinite(bin_width) || bin_width <= 0.0) return;
    for (double value : selected_fwhm) {
        if (!std::isfinite(value)
            || value < diagnostics.histogram_lower
            || value > diagnostics.histogram_upper) {
            continue;
        }
        int bin = static_cast<int>(
            (value - diagnostics.histogram_lower) / bin_width);
        bin = std::clamp(
            bin, 0, static_cast<int>(diagnostics.histogram.size()) - 1);
        diagnostics.selected_group_histogram[static_cast<std::size_t>(bin)]
            += 1.0;
    }
}

// ==========================================
// Function: Parse matched Gaia image positions from one astro stream
// Method: Read the two WCS lines, the matched/user/reference counts, and exactly
//         n_matched finite RA/Dec/x/y rows while retaining only x/y.
// ==========================================
AstrometryGaiaReadStatus parseAstrometryGaiaPositions(
    std::istream& input,
    std::vector<std::array<double, 2>>& gaia_xy,
    std::string& error) {
    gaia_xy.clear();
    error.clear();
    std::string first_line;
    std::string second_line;
    if (!std::getline(input, first_line) || !std::getline(input, second_line)) {
        error = "missing astrometry WCS lines";
        return AstrometryGaiaReadStatus::Malformed;
    }

    int matched = 0;
    int user = 0;
    int reference = 0;
    if (!(input >> matched >> user >> reference) || matched < 0
        || user < 0 || reference < 0) {
        error = "invalid astrometry match-count row";
        return AstrometryGaiaReadStatus::Malformed;
    }

    gaia_xy.reserve(static_cast<std::size_t>(matched));
    for (int row = 0; row < matched; ++row) {
        double ra = 0.0;
        double dec = 0.0;
        double x = 0.0;
        double y = 0.0;
        if (!(input >> ra >> dec >> x >> y)
            || !std::isfinite(ra) || !std::isfinite(dec)
            || !std::isfinite(x) || !std::isfinite(y)) {
            error = "invalid astrometry matched-source row " + std::to_string(row);
            gaia_xy.clear();
            return AstrometryGaiaReadStatus::Malformed;
        }
        gaia_xy.push_back({x, y});
    }
    return gaia_xy.empty()
        ? AstrometryGaiaReadStatus::Empty
        : AstrometryGaiaReadStatus::Accepted;
}

// ==========================================
// Function: Match one candidate to the nearest same-chip Gaia position
// Method: Accept when the minimum squared image-plane distance is within the
//         configured radius; do not impose one-to-one assignment.
// ==========================================
bool hasNearestGaiaMatch(
    double x,
    double y,
    const std::vector<std::array<double, 2>>& gaia_xy,
    double radius_pixels) {
    if (!std::isfinite(x) || !std::isfinite(y)
        || !std::isfinite(radius_pixels) || radius_pixels < 0.0) {
        return false;
    }
    const double radius_squared = radius_pixels * radius_pixels;
    double minimum_squared = std::numeric_limits<double>::infinity();
    for (const std::array<double, 2>& gaia : gaia_xy) {
        const double dx = x - gaia[0];
        const double dy = y - gaia[1];
        minimum_squared = std::min(minimum_squared, dx * dx + dy * dy);
    }
    return minimum_squared <= radius_squared;
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
// Function: Maintain one exact sorted top-K neighbour list
// Method: Insert or improve the candidate edge, sort by chi/index, and truncate.
// ==========================================
void updateTopK(
    std::vector<NeighborEdge>& neighbours,
    int star_index,
    float chi,
    int k) {
    if (star_index < 0 || k <= 0 || !std::isfinite(chi)) return;
    auto existing = std::find_if(
        neighbours.begin(), neighbours.end(),
        [star_index](const NeighborEdge& edge) {
            return edge.star_index == star_index;
        });
    if (existing == neighbours.end()) {
        neighbours.push_back({star_index, chi});
    } else if (chi < existing->chi) {
        existing->chi = chi;
    }
    std::sort(
        neighbours.begin(), neighbours.end(),
        [](const NeighborEdge& first, const NeighborEdge& second) {
            if (first.chi != second.chi) return first.chi < second.chi;
            return first.star_index < second.star_index;
        });
    if (static_cast<int>(neighbours.size()) > k) {
        neighbours.resize(static_cast<std::size_t>(k));
    }
}

// ==========================================
// Function: Select capped exposure-wide large-size minChi references
// Method: Sort the top fraction by size/chip/star, then apply a per-chip cap
//         without filling from candidates below the exposure-wide pool.
// ==========================================
std::vector<MinChiReferenceCandidate> selectMinChiReferenceStars(
    const std::vector<MinChiReferenceCandidate>& locus_candidates,
    double reference_fraction,
    int maximum_per_chip) {
    std::vector<MinChiReferenceCandidate> ranked;
    if (!std::isfinite(reference_fraction) || reference_fraction <= 0.0
        || reference_fraction > 1.0 || maximum_per_chip <= 0) {
        return ranked;
    }
    ranked.reserve(locus_candidates.size());
    for (const MinChiReferenceCandidate& candidate : locus_candidates) {
        if (candidate.chip_index >= 0 && candidate.star_index >= 0
            && std::isfinite(candidate.size)) {
            ranked.push_back(candidate);
        }
    }
    std::sort(
        ranked.begin(), ranked.end(),
        [](const MinChiReferenceCandidate& first,
           const MinChiReferenceCandidate& second) {
            if (first.size != second.size) return first.size > second.size;
            if (first.chip_index != second.chip_index) {
                return first.chip_index < second.chip_index;
            }
            return first.star_index < second.star_index;
        });

    const std::size_t pool_count = std::min(
        ranked.size(),
        static_cast<std::size_t>(std::ceil(
            reference_fraction * static_cast<double>(ranked.size()))));
    std::vector<MinChiReferenceCandidate> selected;
    selected.reserve(pool_count);
    std::vector<int> chip_counts;
    for (std::size_t index = 0; index < pool_count; ++index) {
        const MinChiReferenceCandidate& candidate = ranked[index];
        if (candidate.chip_index >= static_cast<int>(chip_counts.size())) {
            chip_counts.resize(
                static_cast<std::size_t>(candidate.chip_index + 1), 0);
        }
        if (chip_counts[candidate.chip_index] >= maximum_per_chip) continue;
        chip_counts[candidate.chip_index]++;
        selected.push_back(candidate);
    }
    return selected;
}

// ==========================================
// Function: Compute one chip's minChi values and threshold-pair sample
// Method: Visit every unordered locus-locus pair once, update both endpoints,
//         and sample the distance when either endpoint is a reference.
// ==========================================
MinChiPairResult computeMinChiAndThresholdPairs(
    const std::vector<MinChiCandidateView>& candidates) {
    MinChiPairResult result;
    result.min_chi.assign(
        candidates.size(), std::numeric_limits<float>::infinity());
    for (std::size_t first = 0; first + 1 < candidates.size(); ++first) {
        if (!candidates[first].in_size_locus
            || candidates[first].chi_window == nullptr) {
            continue;
        }
        for (std::size_t second = first + 1;
             second < candidates.size(); ++second) {
            if (!candidates[second].in_size_locus
                || candidates[second].chi_window == nullptr) {
                continue;
            }
            const float chi = normalizedChiDistance(
                *candidates[first].chi_window,
                *candidates[second].chi_window);
            if (!std::isfinite(chi)) continue;
            result.min_chi[first] = std::min(result.min_chi[first], chi);
            result.min_chi[second] = std::min(result.min_chi[second], chi);
            if (candidates[first].is_reference
                || candidates[second].is_reference) {
                result.threshold_pair_chi.push_back(chi);
            }
        }
    }
    return result;
}

// ==========================================
// Function: Extract mutual-KNN graph edges among active candidates
// Method: Keep an undirected edge only when both retained top-K lists contain
//         the opposite endpoint and both endpoints survived the shared cut.
// ==========================================
std::vector<GraphEdge> buildMutualKNNEdges(
    const std::vector<int>& active_indices,
    const std::vector<std::vector<NeighborEdge>>& neighbours_by_star) {
    std::vector<GraphEdge> edges;
    if (active_indices.empty()) return edges;
    const int maximum_index = *std::max_element(
        active_indices.begin(), active_indices.end());
    std::vector<bool> active(static_cast<std::size_t>(maximum_index + 1), false);
    for (int index : active_indices) {
        if (index >= 0) active[index] = true;
    }
    for (int first : active_indices) {
        if (first < 0 || first >= static_cast<int>(neighbours_by_star.size())) continue;
        for (const NeighborEdge& neighbour : neighbours_by_star[first]) {
            const int second = neighbour.star_index;
            if (second <= first || second > maximum_index || !active[second]
                || second >= static_cast<int>(neighbours_by_star.size())) {
                continue;
            }
            if (containsNeighbour(neighbours_by_star[second], first)) {
                edges.push_back({first, second});
            }
        }
    }
    return edges;
}

// ==========================================
// Function: Convert one same-chip graph into connected components
// Method: Use disjoint sets over active original indices and attach Gaia counts.
// ==========================================
std::vector<StarGroup> buildConnectedGroups(
    const std::vector<int>& active_indices,
    const std::vector<GraphEdge>& edges,
    const std::vector<bool>& gaia_matched) {
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
        if (original >= 0 && original < static_cast<int>(gaia_matched.size())
            && gaia_matched[original]) {
            group.gaia_count++;
        }
    }
    return groups;
}

// ==========================================
// Function: Select the shared main and eligible secondary stellar groups
// Method: Always keep the largest component and keep another component only
//         when both its relative size and Gaia-count requirements pass.
// ==========================================
std::vector<int> selectMainAndSecondaryGroups(
    const std::vector<StarGroup>& groups,
    double minimum_size_ratio,
    int minimum_gaia_count) {
    std::vector<int> selected;
    if (groups.empty() || minimum_size_ratio < 0.0 || minimum_gaia_count < 0) {
        return selected;
    }
    std::size_t main_group = 0;
    for (std::size_t group = 1; group < groups.size(); ++group) {
        if (groups[group].members.size() > groups[main_group].members.size()) {
            main_group = group;
        }
    }
    const double main_size = static_cast<double>(groups[main_group].members.size());
    for (std::size_t group = 0; group < groups.size(); ++group) {
        const bool keep_main = group == main_group;
        const bool size_ok = main_size > 0.0
            && static_cast<double>(groups[group].members.size())
                >= minimum_size_ratio * main_size;
        const bool gaia_ok = groups[group].gaia_count >= minimum_gaia_count;
        if (keep_main || (size_ok && gaia_ok)) {
            selected.insert(
                selected.end(), groups[group].members.begin(), groups[group].members.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    return selected;
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

// ==========================================
// Function: Convert raw analytic PRESS to its leverage-standardized score
// Method: Multiply by sqrt(1-h) after the same finite denominator guard used
//         by the analytic leave-one-out calculation.
// ==========================================
bool computeLeverageStandardizedPress(
    double raw_press,
    double leverage,
    double minimum_denominator,
    double& standardized_press) {
    const double denominator = 1.0 - leverage;
    if (!std::isfinite(raw_press) || !std::isfinite(denominator)
        || !std::isfinite(minimum_denominator)
        || minimum_denominator <= 0.0
        || denominator <= minimum_denominator) {
        return false;
    }
    standardized_press = raw_press * std::sqrt(denominator);
    return std::isfinite(standardized_press);
}

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
    int maximum_removals) {
    if (!rejection_enabled) return PressRemovalDecision::Disabled;
    if (flagged_count <= 0) return PressRemovalDecision::NoOutliers;
    if (maximum_removals < 0 || flagged_count > maximum_removals) {
        return PressRemovalDecision::TooManyOutliers;
    }
    if (initial_count - flagged_count < minimum_count) {
        return PressRemovalDecision::WouldUnderrunMinimum;
    }
    return PressRemovalDecision::Apply;
}

}  // namespace Internal
}  // namespace PSFModel
