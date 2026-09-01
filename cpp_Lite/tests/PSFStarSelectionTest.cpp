#include "process_main/PSFStarSelection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace PSFModel::Internal;

// ==========================================
// Function: Stop the star-selection test on one failed invariant
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF star-selection test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Stop the test when two finite scalars differ beyond tolerance
// Method: Compare their absolute difference and reuse the focused failure path.
// ==========================================
void requireNear(
    double actual,
    double expected,
    double tolerance,
    const std::string& message) {
    require(std::isfinite(actual)
                && std::abs(actual - expected) <= tolerance,
            message + ": actual=" + std::to_string(actual)
                + ", expected=" + std::to_string(expected));
}

// ==========================================
// Function: Verify shared chi-window bounds and normalized distance
// Method: Check the locked 64-pixel bounds and exact zero/nonzero distances.
// ==========================================
void testChiWindowAndDistance() {
    const PSFChiWindow window = getPSFChiWindow(64);
    require(window.first == 15 && window.last == 47
                && window.pixelCount() == 1089,
            "64x64 chi window must remain inclusive 15..47");
    std::vector<float> first(9, 0.1f);
    std::vector<float> second = first;
    require(normalizedChiDistance(first, second) == 0.0f,
            "identical normalized windows must have zero chi distance");
    second[4] += 0.05f;
    require(normalizedChiDistance(first, second) > 0.0f,
            "different normalized windows must have positive chi distance");
}

// ==========================================
// Function: Build the production-shaped FWHM-locus test configuration
// Method: Keep the production quantile, pilot, and final-cut controls fixed while
//         allowing histogram-resolution checks to vary only the bin count.
// ==========================================
FWHMLocusConfig fwhmLocusConfig(int histogram_bins = 128) {
    return {histogram_bins, 3.0, 3, 0.05, 5.0, 4.0, 30, 10};
}

// ==========================================
// Function: Generate a deterministic smooth single-core FWHM population
// Method: Sample a dense symmetric Gaussian profile on a fine physical grid so
//         histogram-resolution checks are not dominated by random count noise.
// ==========================================
std::vector<FWHMSample> smoothFWHMCore() {
    std::vector<FWHMSample> samples;
    for (int step = -160; step <= 160; ++step) {
        const double offset = 0.00025 * static_cast<double>(step);
        const double normal_density = std::exp(
            -0.5 * (offset / 0.01) * (offset / 0.01));
        const int copies = std::max(
            1, static_cast<int>(std::lround(200.0 * normal_density)));
        for (int copy = 0; copy < copies; ++copy) {
            samples.push_back({1.30 + offset, false});
        }
    }
    return samples;
}

// ==========================================
// Function: Generate a deterministic bin-neutral FWHM population
// Method: Fill one continuous physical interval uniformly so 64/128/256-bin
//         basin boundaries retain the same underlying final population.
// ==========================================
std::vector<FWHMSample> uniformFWHMCore() {
    std::vector<FWHMSample> samples;
    constexpr int sample_count = 2560;
    samples.reserve(sample_count);
    for (int index = 0; index < sample_count; ++index) {
        const double fraction = (static_cast<double>(index) + 0.5)
            / static_cast<double>(sample_count);
        samples.push_back({1.28 + 0.04 * fraction, false});
    }
    return samples;
}

// ==========================================
// Function: Verify same-bin Gaia histogram diagnostics
// Method: Check shape, total accounting, and the per-bin subset invariant.
// ==========================================
void requireGaiaHistogramConsistent(
    const FWHMLocusDiagnostics& diagnostics,
    std::size_t expected_bins,
    const std::string& context) {
    require(diagnostics.gaia_histogram.size() == expected_bins
                && diagnostics.gaia_histogram.size()
                    == diagnostics.histogram.size(),
            context + ": Gaia histogram must use the all-candidate bins");
    require(diagnostics.gaia_histogram_sample_count
                + diagnostics.gaia_histogram_below_count
                + diagnostics.gaia_histogram_above_count
                == diagnostics.gaia_match_count,
            context + ": Gaia histogram counts must cover every Gaia match");
    for (std::size_t bin = 0; bin < diagnostics.histogram.size(); ++bin) {
        require(diagnostics.gaia_histogram[bin] <= diagnostics.histogram[bin],
                context + ": Gaia histogram must be an all-candidate subset");
    }
}

// ==========================================
// Function: Verify quantization-safe first-layer FWHM pilot behavior
// Method: Exercise initial Gaia/all zero MAD, rejected collapsing clips,
//         configurable symmetric quantiles, padding, and Gaia fallback state.
// ==========================================
void testFWHMZeroMadPilot() {
    constexpr double padding = 1.0e-6;
    constexpr double tolerance = 1.0e-12;
    const FWHMLocusConfig config = fwhmLocusConfig();
    FWHMLocus locus;
    FWHMLocusDiagnostics diagnostics;

    std::vector<FWHMSample> gaia_zero_mad;
    gaia_zero_mad.insert(gaia_zero_mad.end(), 60, {2.44, true});
    gaia_zero_mad.insert(gaia_zero_mad.end(), 30, {2.56, true});
    gaia_zero_mad.insert(gaia_zero_mad.end(), 10, {2.70, true});
    require(estimateFWHMLocus(
                gaia_zero_mad, config, locus, &diagnostics),
            "initial zero-MAD Gaia pilot must produce a locus");
    require(diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 100
                && diagnostics.pilot_retained_count == 100
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range
                && !diagnostics.pilot_rejected_zero_mad_clip,
            "initial zero-MAD Gaia pilot must retain its complete population");
    requireNear(diagnostics.pilot_center, 2.44, tolerance,
                "zero-MAD Gaia center must remain the initial median");
    requireNear(diagnostics.pilot_lower, 2.44 - padding, tolerance,
                "zero-MAD Gaia lower bound must be padded Q05");
    requireNear(diagnostics.pilot_upper, 2.70 + padding, tolerance,
                "zero-MAD Gaia upper bound must be padded Q95");
    requireGaiaHistogramConsistent(
        diagnostics, 128, "zero-MAD Gaia fixture");

    std::vector<FWHMSample> all_zero_mad;
    all_zero_mad.insert(all_zero_mad.end(), 70, {1.25, false});
    all_zero_mad.insert(all_zero_mad.end(), 20, {1.35, false});
    all_zero_mad.insert(all_zero_mad.end(), 10, {1.50, false});
    require(estimateFWHMLocus(
                all_zero_mad, config, locus, &diagnostics),
            "initial zero-MAD all-candidate pilot must produce a locus");
    require(!diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 100
                && diagnostics.pilot_retained_count == 100
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range
                && !diagnostics.pilot_rejected_zero_mad_clip,
            "initial zero-MAD all-candidate pilot must skip clipping");
    requireNear(diagnostics.pilot_lower, 1.25 - padding, tolerance,
                "zero-MAD all-candidate lower bound must be padded Q05");
    requireNear(diagnostics.pilot_upper, 1.50 + padding, tolerance,
                "zero-MAD all-candidate upper bound must be padded Q95");

    FWHMLocusConfig custom_quantile_config = config;
    custom_quantile_config.zero_mad_quantile = 0.20;
    std::vector<FWHMSample> custom_quantile_samples;
    for (int index = 0; index < 24; ++index) {
        custom_quantile_samples.push_back({0.76 + 0.01 * index, false});
    }
    custom_quantile_samples.insert(
        custom_quantile_samples.end(), 52, {1.00, false});
    for (int index = 1; index <= 24; ++index) {
        custom_quantile_samples.push_back({1.00 + 0.01 * index, false});
    }
    require(estimateFWHMLocus(
                custom_quantile_samples,
                custom_quantile_config,
                locus,
                &diagnostics),
            "a configurable Q20--Q80 zero-MAD pilot must produce a locus");
    require(diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range,
            "the custom quantile fixture must exercise the initial zero-MAD branch");
    requireNear(diagnostics.pilot_lower, 0.958 - padding, tolerance,
                "custom zero-MAD lower bound must use interpolated Q20");
    requireNear(diagnostics.pilot_upper, 1.042 + padding, tolerance,
                "custom zero-MAD upper bound must use interpolated Q80");

    FWHMLocusConfig invalid_quantile_config = config;
    invalid_quantile_config.zero_mad_quantile = -0.01;
    require(!estimateFWHMLocus(
                custom_quantile_samples, invalid_quantile_config, locus),
            "a negative zero-MAD quantile must be rejected");
    invalid_quantile_config.zero_mad_quantile = 0.50;
    require(!estimateFWHMLocus(
                custom_quantile_samples, invalid_quantile_config, locus),
            "a half-or-greater zero-MAD quantile must be rejected");

    std::vector<FWHMSample> collapsing_clip;
    collapsing_clip.insert(collapsing_clip.end(), 20, {0.90, false});
    collapsing_clip.insert(collapsing_clip.end(), 45, {1.00, false});
    collapsing_clip.insert(collapsing_clip.end(), 20, {1.10, false});
    collapsing_clip.insert(collapsing_clip.end(), 15, {3.00, false});
    require(estimateFWHMLocus(
                collapsing_clip, config, locus, &diagnostics),
            "a proposed zero-MAD clipping step must be recoverable");
    const double expected_width = 1.4826 * 0.10;
    require(!diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 100
                && diagnostics.pilot_retained_count == 100
                && diagnostics.pilot_width > 0.0
                && !diagnostics.pilot_uses_quantile_range
                && diagnostics.pilot_rejected_zero_mad_clip,
            "a collapsing clip must be rejected before population swap");
    requireNear(diagnostics.pilot_center, 1.00, tolerance,
                "rejected clip must retain the previous pilot center");
    requireNear(diagnostics.pilot_width, expected_width, tolerance,
                "rejected clip must retain the previous positive MAD");
    requireNear(
        diagnostics.pilot_lower,
        1.00 - config.histogram_range_sigma * expected_width - padding,
        tolerance,
        "normal pilot lower bound must include outward padding");
    requireNear(
        diagnostics.pilot_upper,
        1.00 + config.histogram_range_sigma * expected_width + padding,
        tolerance,
        "normal pilot upper bound must include outward padding");

    const FWHMLocus default_positive_locus = locus;
    const FWHMLocusDiagnostics default_positive_diagnostics = diagnostics;
    require(estimateFWHMLocus(
                collapsing_clip,
                custom_quantile_config,
                locus,
                &diagnostics),
            "the custom quantile must leave a positive-MAD pilot valid");
    require(!diagnostics.pilot_uses_quantile_range
                && diagnostics.pilot_center
                    == default_positive_diagnostics.pilot_center
                && diagnostics.pilot_width
                    == default_positive_diagnostics.pilot_width
                && diagnostics.pilot_lower
                    == default_positive_diagnostics.pilot_lower
                && diagnostics.pilot_upper
                    == default_positive_diagnostics.pilot_upper
                && locus.center == default_positive_locus.center
                && locus.width == default_positive_locus.width
                && locus.lower == default_positive_locus.lower
                && locus.upper == default_positive_locus.upper,
            "the zero-MAD quantile must not affect positive-MAD behavior");

    std::vector<FWHMSample> fallback_zero_mad;
    for (int index = 0; index < 9; ++index) {
        fallback_zero_mad.push_back({1.30 + 0.002 * index, true});
    }
    fallback_zero_mad.push_back({2.50, true});
    fallback_zero_mad.insert(fallback_zero_mad.end(), 70, {1.25, false});
    fallback_zero_mad.insert(fallback_zero_mad.end(), 20, {1.35, false});
    require(estimateFWHMLocus(
                fallback_zero_mad, config, locus, &diagnostics),
            "Gaia retained-count fallback to a zero-MAD all pilot must succeed");
    require(diagnostics.has_gaia_median
                && !diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 100
                && diagnostics.pilot_retained_count == 100
                && diagnostics.pilot_center == 1.25
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range
                && !diagnostics.pilot_rejected_zero_mad_clip,
            "Gaia fallback must replace every pilot field and branch flag");
    requireNear(diagnostics.pilot_lower, 1.25 - padding, tolerance,
                "fallback pilot lower bound must come from all-candidate Q05");
    requireNear(diagnostics.pilot_upper, 1.35 + padding, tolerance,
                "fallback pilot upper bound must come from all-candidate Q95");
}

// ==========================================
// Function: Verify robust-pilot FWHM peak selection and width independence
// Method: Exercise density/Gaia routes, long tails, Gaia fallback, repeated
//         values, window accounting, bin invariance, and sample-count failure.
// ==========================================
void testFWHMLocus() {
    constexpr double padding = 1.0e-6;
    constexpr double tolerance = 1.0e-12;
    const FWHMLocusConfig standard_config = fwhmLocusConfig();
    std::vector<FWHMSample> density_samples;
    for (int index = 0; index < 80; ++index) {
        density_samples.push_back({1.00 + 0.01 * (index % 5), false});
    }
    for (int index = 0; index < 120; ++index) {
        density_samples.push_back({1.70 + 0.02 * (index % 25), false});
    }
    FWHMLocus locus;
    FWHMLocusDiagnostics diagnostics;
    require(estimateFWHMLocus(
                density_samples, standard_config, locus, &diagnostics),
            "two-population FWHM fixture must produce a locus");
    require(std::abs(locus.center - 1.02) < 0.08,
            "highest-density narrow stellar peak must be selected without Gaia");
    require(diagnostics.sample_count == 200
                && diagnostics.histogram.size() == 128
                && diagnostics.smoothed_histogram.size() == 128
                && diagnostics.peak_bin >= 0 && diagnostics.peak_bin < 128
                && diagnostics.pilot_upper > diagnostics.pilot_lower
                && diagnostics.histogram_sample_count
                    + diagnostics.histogram_below_count
                    + diagnostics.histogram_above_count
                    == diagnostics.sample_count,
            "density fixture must publish its exact histogram diagnostics");
    require(!diagnostics.pilot_uses_gaia && !diagnostics.has_gaia_median,
            "density-only peak selection must use the all-candidate pilot");
    requireGaiaHistogramConsistent(diagnostics, 128, "zero-Gaia fixture");
    require(diagnostics.gaia_match_count == 0
                && diagnostics.gaia_histogram_sample_count == 0
                && std::all_of(
                    diagnostics.gaia_histogram.begin(),
                    diagnostics.gaia_histogram.end(),
                    [](double count) { return count == 0.0; }),
            "zero-Gaia fixture must publish a stable all-zero histogram");
    FWHMLocus baseline_locus;
    require(estimateFWHMLocus(
                density_samples, standard_config, baseline_locus)
                && baseline_locus.center == locus.center
                && baseline_locus.width == locus.width
                && baseline_locus.lower == locus.lower
                && baseline_locus.upper == locus.upper,
            "requesting diagnostics must not alter the scientific locus");

    std::vector<FWHMSample> gaia_samples;
    for (int index = 0; index < 45; ++index) {
        gaia_samples.push_back({0.90 + 0.01 * (index % 5), index < 12});
    }
    for (int index = 0; index < 100; ++index) {
        gaia_samples.push_back({1.80 + 0.005 * (index % 5), false});
    }
    require(estimateFWHMLocus(
                gaia_samples, standard_config, locus, &diagnostics),
            "Gaia-supported two-peak fixture must produce a locus");
    require(std::abs(locus.center - 0.92) < 0.08,
            "Gaia median must select the supported smaller peak");
    require(diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 12
                && diagnostics.pilot_retained_count == 12
                && diagnostics.has_gaia_median
                && diagnostics.gaia_match_count == 12
                && std::isfinite(diagnostics.gaia_median),
            "Gaia-supported selection must publish raw and clipped-pilot diagnostics");
    requireGaiaHistogramConsistent(diagnostics, 128, "Gaia-supported fixture");

    std::vector<FWHMSample> gaia_tail_samples;
    for (int index = 0; index < 100; ++index) {
        gaia_tail_samples.push_back({
            1.28 + 0.002 * static_cast<double>(index % 21), index < 20});
    }
    for (int index = 0; index < 20; ++index) {
        gaia_tail_samples.push_back({2.0 + 0.15 * index, false});
    }
    require(estimateFWHMLocus(
                gaia_tail_samples, standard_config, locus, &diagnostics),
            "Gaia-supported long-tail fixture must produce a locus");
    require(diagnostics.pilot_uses_gaia
                && diagnostics.histogram_above_count > 0
                && diagnostics.pilot_upper - diagnostics.pilot_lower < 0.5
                && locus.upper < 1.5,
            "Gaia pilot must keep the high-FWHM tail out of the local window");

    std::vector<FWHMSample> all_tail_samples;
    for (int index = 0; index < 140; ++index) {
        all_tail_samples.push_back({
            1.25 + 0.002 * static_cast<double>(index % 25), index < 5});
    }
    for (int index = 0; index < 15; ++index) {
        all_tail_samples.push_back({2.0 + 0.2 * index, false});
    }
    require(estimateFWHMLocus(
                all_tail_samples, standard_config, locus, &diagnostics),
            "all-candidate long-tail fixture must produce a locus");
    require(!diagnostics.pilot_uses_gaia
                && diagnostics.gaia_match_count == 5
                && diagnostics.pilot_input_count == 155
                && diagnostics.pilot_retained_count < 155
                && diagnostics.histogram_above_count > 0
                && diagnostics.pilot_upper - diagnostics.pilot_lower < 0.5,
            "all-candidate pilot must iteratively clip and exclude its high tail");
    requireGaiaHistogramConsistent(diagnostics, 128, "insufficient-Gaia fixture");
    require(diagnostics.gaia_histogram_sample_count > 0,
            "insufficient Gaia must still publish its local-window distribution");

    std::vector<FWHMSample> gaia_clipped_samples;
    for (int index = 0; index < 10; ++index) {
        gaia_clipped_samples.push_back({1.30 + 0.002 * index, true});
    }
    gaia_clipped_samples.push_back({2.50, true});
    for (int index = 0; index < 60; ++index) {
        gaia_clipped_samples.push_back({
            1.28 + 0.001 * static_cast<double>(index % 41), false});
    }
    require(estimateFWHMLocus(
                gaia_clipped_samples, standard_config, locus, &diagnostics),
            "retained-sufficient clipped Gaia fixture must produce a locus");
    require(diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 11
                && diagnostics.pilot_retained_count == 10,
            "Gaia pilot must apply the same 3-MAD clipping while retaining support");

    std::vector<FWHMSample> gaia_fallback_samples;
    for (int index = 0; index < 9; ++index) {
        gaia_fallback_samples.push_back({1.30 + 0.002 * index, true});
    }
    gaia_fallback_samples.push_back({2.50, true});
    for (int index = 0; index < 90; ++index) {
        gaia_fallback_samples.push_back({
            1.28 + 0.001 * static_cast<double>(index % 41), false});
    }
    require(estimateFWHMLocus(
                gaia_fallback_samples, standard_config, locus, &diagnostics),
            "clipped-insufficient Gaia fixture must fall back and produce a locus");
    require(diagnostics.has_gaia_median && !diagnostics.pilot_uses_gaia
                && diagnostics.pilot_input_count == 100
                && diagnostics.pilot_retained_count < 100,
            "Gaia pilot below its retained minimum must rerun on all candidates");

    std::vector<FWHMSample> repeated(40, {1.25, false});
    for (int index = 0; index < 6; ++index) {
        repeated[static_cast<std::size_t>(index)].gaia_matched = true;
    }
    require(estimateFWHMLocus(
                repeated, standard_config, locus, &diagnostics)
                && locus.width > 0.0 && locus.lower < 1.25 && locus.upper > 1.25,
            "repeated FWHM values must receive a finite positive width floor");
    require(diagnostics.histogram.size() == 128
                && diagnostics.smoothed_histogram.size() == 128
                && diagnostics.peak_bin >= 0
                && diagnostics.pilot_width == 0.0
                && diagnostics.pilot_uses_quantile_range
                && !diagnostics.pilot_rejected_zero_mad_clip
                && diagnostics.pilot_upper > diagnostics.pilot_lower
                && diagnostics.histogram_sample_count == 40
                && diagnostics.histogram_below_count == 0
                && diagnostics.histogram_above_count == 0,
            "Q05==Q95 must use padded bounds without a special one-bin branch");
    requireNear(diagnostics.pilot_lower, 1.25 - padding, tolerance,
                "repeated-value lower bound must be padded Q05");
    requireNear(diagnostics.pilot_upper, 1.25 + padding, tolerance,
                "repeated-value upper bound must be padded Q95");
    requireGaiaHistogramConsistent(diagnostics, 128, "repeated-value fixture");
    require(diagnostics.gaia_histogram_sample_count == 6
                && std::count(
                    diagnostics.gaia_histogram.begin(),
                    diagnostics.gaia_histogram.end(), 6.0) == 1,
            "repeated values must retain all Gaia matches in one shared bin");

    const std::vector<FWHMSample> invariant_samples = uniformFWHMCore();
    FWHMLocus locus_64;
    FWHMLocus locus_128;
    FWHMLocus locus_256;
    require(estimateFWHMLocus(
                invariant_samples, fwhmLocusConfig(64), locus_64)
                && estimateFWHMLocus(
                    invariant_samples, fwhmLocusConfig(128), locus_128)
                && estimateFWHMLocus(
                    invariant_samples, fwhmLocusConfig(256), locus_256),
            "bin-neutral core must produce loci for 64, 128, and 256 bins");
    require(std::abs(locus_64.center - locus_128.center) < 0.002
                && std::abs(locus_256.center - locus_128.center) < 0.002
                && std::abs(locus_64.width - locus_128.width) < 0.002
                && std::abs(locus_256.width - locus_128.width) < 0.002
                && std::abs(locus_64.lower - locus_128.lower) < 0.01
                && std::abs(locus_256.lower - locus_128.lower) < 0.01
                && std::abs(locus_64.upper - locus_128.upper) < 0.01
                && std::abs(locus_256.upper - locus_128.upper) < 0.01,
            "final locus must remain stable across histogram bin counts: "
                + std::to_string(locus_64.center) + ","
                + std::to_string(locus_64.width) + " / "
                + std::to_string(locus_128.center) + ","
                + std::to_string(locus_128.width) + " / "
                + std::to_string(locus_256.center) + ","
                + std::to_string(locus_256.width));
    const std::vector<FWHMSample> smooth_samples = smoothFWHMCore();
    FWHMLocus coarse_locus;
    require(estimateFWHMLocus(
                smooth_samples, fwhmLocusConfig(3), coarse_locus)
                && coarse_locus.width < coarse_locus.histogram_bin_width,
            "final physical width must be allowed below a coarse histogram bin");

    repeated.resize(29);
    require(!estimateFWHMLocus(
                repeated, standard_config, locus, &diagnostics),
            "FWHM locus must reject fewer than the configured samples");
    require(diagnostics.sample_count == 29 && diagnostics.histogram.empty(),
            "failed estimation must reset diagnostics before reporting counts");
}

// ==========================================
// Function: Verify astro parsing and nearest-position Gaia labels
// Method: Parse finite matched rows, test radius boundaries and shared labels,
//         and reject a malformed matched-source row.
// ==========================================
void testGaiaParsingAndMatching() {
    std::istringstream valid(
        "1 2 3 4\n"
        "0.1 0.0 0.0 0.1\n"
        "2 20 30\n"
        "10 20 100.0 200.0\n"
        "11 21 150.0 250.0\n");
    std::vector<std::array<double, 2>> gaia_xy;
    std::string error;
    require(parseAstrometryGaiaPositions(valid, gaia_xy, error)
                == AstrometryGaiaReadStatus::Accepted,
            "well-formed astro rows must parse");
    require(gaia_xy.size() == 2
                && hasNearestGaiaMatch(101.0, 201.0, gaia_xy, 2.5)
                && !hasNearestGaiaMatch(104.0, 204.0, gaia_xy, 2.5),
            "nearest Gaia radius must be applied in image pixels");
    require(hasNearestGaiaMatch(100.5, 200.5, gaia_xy, 2.5)
                && hasNearestGaiaMatch(99.5, 199.5, gaia_xy, 2.5),
            "multiple candidates may share one supporting Gaia position");

    std::istringstream malformed(
        "1 2 3 4\n0.1 0.0 0.0 0.1\n1 2 3\n10 20 bad 40\n");
    require(parseAstrometryGaiaPositions(malformed, gaia_xy, error)
                == AstrometryGaiaReadStatus::Malformed,
            "malformed astro matched-source rows must be rejected");
}

// ==========================================
// Function: Verify streaming top-K, mutual edges, and shared group selection
// Method: Compare deterministic top-K order, build mutual components, and
//         exercise the secondary size-and-Gaia conjunction.
// ==========================================
void testGrouping() {
    std::vector<NeighborEdge> top_k;
    updateTopK(top_k, 3, 3.0f, 2);
    updateTopK(top_k, 1, 1.0f, 2);
    updateTopK(top_k, 2, 2.0f, 2);
    require(top_k.size() == 2 && top_k[0].star_index == 1
                && top_k[1].star_index == 2,
            "streaming top-K must match sorted full-distance reference");

    std::vector<std::vector<NeighborEdge>> neighbours(4);
    neighbours[0] = {{1, 1.0f}, {2, 2.0f}};
    neighbours[1] = {{0, 1.0f}, {2, 1.5f}};
    neighbours[2] = {{1, 1.5f}, {0, 2.0f}};
    neighbours[3] = {{2, 0.5f}};
    const std::vector<int> active = {0, 1, 2, 3};
    const std::vector<GraphEdge> mutual = buildMutualKNNEdges(active, neighbours);
    const std::vector<bool> gaia = {false, true, false, true};
    const std::vector<StarGroup> connected =
        buildConnectedGroups(active, mutual, gaia);
    require(connected.size() == 2,
            "mutual-KNN fixture must form one triple and one singleton");

    const std::vector<StarGroup> groups = {
        {{0, 1, 2, 3}, 0},
        {{4, 5}, 2},
        {{6, 7}, 1},
        {{8}, 3}
    };
    const std::vector<int> selected =
        selectMainAndSecondaryGroups(groups, 0.40, 2);
    require(selected == std::vector<int>({0, 1, 2, 3, 4, 5}),
            "secondary groups must pass both relative size and Gaia count");
}

// ==========================================
// Structure: Provide the minimal cached fields required by active-KNN rebuilds
// Method: Mirror production chi-window and neighbour storage without PSF I/O.
// ==========================================
struct SyntheticKNNSelectionState {
    std::vector<float> chi_window;
    std::vector<NeighborEdge> knn;
};

// ==========================================
// Function: Verify KNN slots are refilled after the minChi survivor cut
// Method: Build an initial top-2 containing rejected close stars, rebuild on
//         survivors only, and require the next two valid neighbours to replace them.
// ==========================================
void testKNNRebuiltAfterMinChiCut() {
    std::vector<SyntheticKNNSelectionState> candidates(5);
    for (int index = 0; index < static_cast<int>(candidates.size()); ++index) {
        candidates[index].chi_window = {
            1.0f + static_cast<float>(index) * 0.1f};
    }

    rebuildActiveKNN(std::vector<int>({0, 1, 2, 3, 4}), candidates, 2);
    require(candidates[0].knn.size() == 2
                && candidates[0].knn[0].star_index == 1
                && candidates[0].knn[1].star_index == 2,
            "pre-cut top-K fixture must be occupied by the closest rejected stars");

    rebuildActiveKNN(std::vector<int>({0, 3, 4}), candidates, 2);
    require(candidates[0].knn.size() == 2
                && candidates[0].knn[0].star_index == 3
                && candidates[0].knn[1].star_index == 4,
            "survivor-only rebuild must refill every vacated top-K slot");
    require(candidates[1].knn.empty() && candidates[2].knn.empty(),
            "minChi-rejected candidates must retain no stale neighbour state");
}

// ==========================================
// Function: Verify capped reference selection and one-pass minChi semantics
// Method: Exercise exposure top-fraction ranking, per-chip caps, deterministic
//         ties, non-reference nearest pairs, and unordered-pair deduplication.
// ==========================================
void testMinChiReferencesAndPairs() {
    std::vector<MinChiReferenceCandidate> candidates;
    for (int star = 0; star < 6; ++star) {
        candidates.push_back({0, star, 100.0 - star});
    }
    candidates.push_back({1, 0, 95.0});
    for (int star = 1; star < 15; ++star) {
        candidates.push_back({1, star, 94.0 - star});
    }
    const std::vector<MinChiReferenceCandidate> selected =
        selectMinChiReferenceStars(candidates, 1.0 / 3.0, 5);
    require(selected.size() == 6,
            "top-third reference pool must apply the cap without backfilling");
    for (int star = 0; star < 5; ++star) {
        require(selected[star].chip_index == 0
                    && selected[star].star_index == star,
                "largest five candidates on one chip must be retained in order");
    }
    require(selected[5].chip_index == 1 && selected[5].star_index == 0,
            "a capped top-pool entry must not suppress another chip reference");

    const std::vector<MinChiReferenceCandidate> tied = {
        {1, 2, 10.0}, {0, 5, 10.0}, {0, 3, 10.0}, {1, 1, 10.0},
        {2, 0, 9.0}, {2, 1, 8.0}, {2, 2, 7.0}, {2, 3, 6.0}, {2, 4, 5.0}
    };
    const std::vector<MinChiReferenceCandidate> tied_selected =
        selectMinChiReferenceStars(tied, 1.0 / 3.0, 5);
    require(tied_selected.size() == 3
                && tied_selected[0].chip_index == 0
                && tied_selected[0].star_index == 3
                && tied_selected[1].chip_index == 0
                && tied_selected[1].star_index == 5
                && tied_selected[2].chip_index == 1
                && tied_selected[2].star_index == 1,
            "size ties must resolve by chip index and then star index");

    const std::vector<std::vector<float>> windows = {
        {1.0f}, {1.2f}, {5.0f}, {5.01f}
    };
    std::vector<MinChiCandidateView> views;
    for (int index = 0; index < 4; ++index) {
        views.push_back({&windows[index], true, index < 2});
    }
    const MinChiPairResult pair_result =
        computeMinChiAndThresholdPairs(views);
    const float nonreference_distance =
        normalizedChiDistance(windows[2], windows[3]);
    require(std::abs(pair_result.min_chi[2] - nonreference_distance) < 1.0e-7f
                && std::abs(pair_result.min_chi[3] - nonreference_distance)
                    < 1.0e-7f,
            "a non-reference pair must still update both nearest distances");
    require(pair_result.threshold_pair_chi.size() == 5,
            "two references among four stars must contribute five unique pairs");
    const float reference_pair_distance =
        normalizedChiDistance(windows[0], windows[1]);
    int reference_pair_occurrences = 0;
    for (float chi : pair_result.threshold_pair_chi) {
        if (std::abs(chi - reference_pair_distance) < 1.0e-7f) {
            reference_pair_occurrences++;
        }
    }
    require(reference_pair_occurrences == 1,
            "the reference-reference pair must enter the threshold sample once");
}

// ==========================================
// Function: Verify leverage-standardized PRESS and removal safeguards
// Method: Check the exact score formula, denominator guard, removal cap, and
//         inclusive retained-star minimum boundaries.
// ==========================================
void testPressStandardizationAndDecision() {
    double standardized = 0.0;
    require(computeLeverageStandardizedPress(
                4.0, 0.75, 1.0e-6, standardized)
                && std::abs(standardized - 2.0) < 1.0e-12,
            "standardized PRESS must equal raw PRESS times sqrt(1-h)");
    require(!computeLeverageStandardizedPress(
                4.0, 1.0 - 5.0e-7, 1.0e-6, standardized),
            "standardized PRESS must reject a denominator below its floor");
    require(!computeLeverageStandardizedPress(
                std::numeric_limits<double>::infinity(), 0.1, 1.0e-6,
                standardized),
            "standardized PRESS must reject non-finite raw scores");

    require(decidePressRemoval(false, 30, 2, 16, 5)
                == PressRemovalDecision::Disabled,
            "disabled rejection must preserve the first fit");
    require(decidePressRemoval(true, 30, 0, 16, 5)
                == PressRemovalDecision::NoOutliers,
            "zero outliers must preserve the first fit");
    require(decidePressRemoval(true, 21, 5, 16, 5)
                == PressRemovalDecision::Apply,
            "five removals leaving sixteen stars must be allowed");
    require(decidePressRemoval(true, 30, 6, 16, 5)
                == PressRemovalDecision::TooManyOutliers,
            "six removals must trigger the configured cap");
    require(decidePressRemoval(true, 20, 5, 16, 5)
                == PressRemovalDecision::WouldUnderrunMinimum,
            "fifteen retained stars must trigger the local-minimum guard");
}

// ==========================================
// Function: Invert one finite nonsingular three-by-three matrix
// Method: Apply pivoted Gauss-Jordan elimination to a small augmented matrix.
// ==========================================
bool invertThreeByThree(
    const std::array<std::array<double, 3>, 3>& matrix,
    std::array<std::array<double, 3>, 3>& inverse) {
    std::array<std::array<double, 6>, 3> augmented{};
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            augmented[row][column] = matrix[row][column];
        }
        augmented[row][row + 3] = 1.0;
    }
    for (int pivot = 0; pivot < 3; ++pivot) {
        int best_row = pivot;
        for (int row = pivot + 1; row < 3; ++row) {
            if (std::abs(augmented[row][pivot])
                > std::abs(augmented[best_row][pivot])) {
                best_row = row;
            }
        }
        if (std::abs(augmented[best_row][pivot]) < 1.0e-14) return false;
        std::swap(augmented[pivot], augmented[best_row]);
        const double scale = augmented[pivot][pivot];
        for (double& value : augmented[pivot]) value /= scale;
        for (int row = 0; row < 3; ++row) {
            if (row == pivot) continue;
            const double factor = augmented[row][pivot];
            for (int column = 0; column < 6; ++column) {
                augmented[row][column] -= factor * augmented[pivot][column];
            }
        }
    }
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            inverse[row][column] = augmented[row][column + 3];
        }
    }
    return true;
}

// ==========================================
// Function: Fit a quadratic while optionally omitting one sample
// Method: Form the three-term normal equations, invert them, and multiply by
//         the observation cross-products for the dependency-free LOO fixture.
// ==========================================
bool fitQuadratic(
    const std::vector<double>& x,
    const std::vector<double>& observed,
    int omitted,
    std::array<double, 3>& coefficients,
    std::array<std::array<double, 3>, 3>* covariance = nullptr) {
    std::array<std::array<double, 3>, 3> normal{};
    std::array<double, 3> rhs{};
    for (int row = 0; row < static_cast<int>(x.size()); ++row) {
        if (row == omitted) continue;
        const std::array<double, 3> basis = {1.0, x[row], x[row] * x[row]};
        for (int first = 0; first < 3; ++first) {
            rhs[first] += basis[first] * observed[row];
            for (int second = 0; second < 3; ++second) {
                normal[first][second] += basis[first] * basis[second];
            }
        }
    }
    std::array<std::array<double, 3>, 3> inverse{};
    if (!invertThreeByThree(normal, inverse)) return false;
    for (int first = 0; first < 3; ++first) {
        for (int second = 0; second < 3; ++second) {
            coefficients[first] += inverse[first][second] * rhs[second];
        }
    }
    if (covariance != nullptr) *covariance = inverse;
    return true;
}

// ==========================================
// Function: Verify analytic LOO against explicit leave-one-out fits
// Method: Fit a small quadratic design, use each hat diagonal in the production
//         formula, and compare with a brute-force refit that removes that row.
// ==========================================
void testAnalyticLOO() {
    constexpr int sample_count = 20;
    std::vector<double> x(static_cast<std::size_t>(sample_count));
    std::vector<double> observed(static_cast<std::size_t>(sample_count));
    for (int row = 0; row < sample_count; ++row) {
        x[row] = -1.0 + 2.0 * row / static_cast<double>(sample_count - 1);
        observed[row] = 2.0 - 0.4 * x[row] + 0.7 * x[row] * x[row]
            + 0.02 * std::sin(3.0 * row);
    }
    std::array<double, 3> coefficients{};
    std::array<std::array<double, 3>, 3> covariance{};
    require(fitQuadratic(x, observed, -1, coefficients, &covariance),
            "full quadratic fixture must be nonsingular");

    for (int omitted = 0; omitted < sample_count; ++omitted) {
        const std::array<double, 3> basis = {
            1.0, x[omitted], x[omitted] * x[omitted]};
        double fitted = 0.0;
        double leverage = 0.0;
        for (int first = 0; first < 3; ++first) {
            fitted += basis[first] * coefficients[first];
            for (int second = 0; second < 3; ++second) {
                leverage += basis[first] * covariance[first][second]
                    * basis[second];
            }
        }
        double loo_residual = 0.0;
        double loo_model = 0.0;
        require(computeAnalyticLOO(
                    observed[omitted], fitted, leverage,
                    1.0e-6, loo_residual, loo_model),
                "well-conditioned leverage must produce an analytic LOO model");

        std::array<double, 3> reduced_coefficients{};
        require(fitQuadratic(
                    x, observed, omitted, reduced_coefficients, nullptr),
                "leave-one-out quadratic fixture must remain nonsingular");
        double explicit_model = 0.0;
        for (int term = 0; term < 3; ++term) {
            explicit_model += basis[term] * reduced_coefficients[term];
        }
        require(std::abs(loo_model - explicit_model) < 1.0e-10,
                "analytic and explicit leave-one-out predictions must agree");
        require(std::abs(loo_residual - (observed[omitted] - explicit_model))
                    < 1.0e-10,
                "analytic LOO residual must equal observed minus explicit model");
    }

    double residual = 0.0;
    double model = 0.0;
    require(!computeAnalyticLOO(1.0, 0.0, 1.0 - 5.0e-7,
                                1.0e-6, residual, model),
            "near-unit leverage must fail the configured denominator guard");
}

}  // namespace

// ==========================================
// Function: Run the focused PSF star-selection regression suite
// Method: Execute quality-independent locus, Gaia, grouping, and LOO cases.
// ==========================================
int main() {
    testChiWindowAndDistance();
    testFWHMZeroMadPilot();
    testFWHMLocus();
    testGaiaParsingAndMatching();
    testGrouping();
    testKNNRebuiltAfterMinChiCut();
    testMinChiReferencesAndPairs();
    testAnalyticLOO();
    testPressStandardizationAndDecision();
    std::cout << "PSF star-selection tests passed\n";
    return EXIT_SUCCESS;
}
