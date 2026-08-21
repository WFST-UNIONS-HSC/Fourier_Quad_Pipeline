#include "NoiseCovariance.hpp"
#include "LensingConfig.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

    // ==========================================
    // Function: Map a signed lag into a compact test covariance grid
    // Method: Mirror the production (2*L+1)^2 row-major lag layout.
    // ==========================================
    std::size_t lagIndex(int dx, int dy, int maxLag) {
        const int side = 2 * maxLag + 1;
        return static_cast<std::size_t>(dy + maxLag) * static_cast<std::size_t>(side)
             + static_cast<std::size_t>(dx + maxLag);
    }

    // ==========================================
    // Function: Compute the masked covariance by explicit pixel-pair summation
    // Method: Enumerate every non-wrapping pair at each signed lag and apply the same minimum
    //         pair-fraction and reversal symmetry contract as the FFT implementation.
    // ==========================================
    std::vector<double> directCovariance(
        int regionSize, int maxLag, double minPairFraction,
        const std::vector<double>& residual,
        const std::vector<unsigned char>& mask) {
        const int side = 2 * maxLag + 1;
        std::vector<double> covariance(static_cast<std::size_t>(side) * side, 0.0);
        std::vector<double> pairCounts(static_cast<std::size_t>(side) * side, 0.0);
        double zeroLagPairs = 0.0;
        for (unsigned char value : mask) {
            if (value != 0U) zeroLagPairs += 1.0;
        }
        const double minimumPairs = minPairFraction * zeroLagPairs;

        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                double numerator = 0.0;
                double pairs = 0.0;
                for (int y = 0; y < regionSize; ++y) {
                    const int otherY = y + dy;
                    if (otherY < 0 || otherY >= regionSize) continue;
                    for (int x = 0; x < regionSize; ++x) {
                        const int otherX = x + dx;
                        if (otherX < 0 || otherX >= regionSize) continue;
                        const std::size_t first = static_cast<std::size_t>(y) * regionSize + x;
                        const std::size_t second = static_cast<std::size_t>(otherY)
                                                 * regionSize + otherX;
                        if (mask[first] == 0U || mask[second] == 0U) continue;
                        numerator += residual[first] * residual[second];
                        pairs += 1.0;
                    }
                }
                const std::size_t index = lagIndex(dx, dy, maxLag);
                pairCounts[index] = pairs;
                if (pairs + 1.0e-8 >= minimumPairs) {
                    covariance[index] = numerator / pairs;
                }
            }
        }

        for (int dy = 0; dy <= maxLag; ++dy) {
            const int dxStart = dy == 0 ? 0 : -maxLag;
            for (int dx = dxStart; dx <= maxLag; ++dx) {
                const std::size_t forward = lagIndex(dx, dy, maxLag);
                const std::size_t reverse = lagIndex(-dx, -dy, maxLag);
                if (pairCounts[forward] + 1.0e-8 >= minimumPairs
                    && pairCounts[reverse] + 1.0e-8 >= minimumPairs) {
                    const double average = 0.5 * (covariance[forward] + covariance[reverse]);
                    covariance[forward] = average;
                    covariance[reverse] = average;
                } else {
                    covariance[forward] = 0.0;
                    covariance[reverse] = 0.0;
                }
            }
        }
        return covariance;
    }

    // ==========================================
    // Function: Compare two numeric vectors within scaled floating-point tolerance
    // Method: Use a mixed absolute/relative error bound and report the first mismatch.
    // ==========================================
    bool vectorsNear(const std::vector<double>& first,
                     const std::vector<double>& second,
                     double tolerance,
                     const char* label) {
        if (first.size() != second.size()) {
            std::cerr << label << ": size mismatch\n";
            return false;
        }
        for (std::size_t i = 0; i < first.size(); ++i) {
            const double scale = std::max({1.0, std::abs(first[i]), std::abs(second[i])});
            if (std::abs(first[i] - second[i]) > tolerance * scale) {
                std::cerr << label << ": mismatch at " << i << " ("
                          << first[i] << " vs " << second[i] << ")\n";
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Compare FFT and direct covariance for one mask geometry
    // Method: Run both estimators on identical deterministic residuals and require agreement.
    // ==========================================
    bool compareCovariance(const std::vector<double>& residual,
                           const std::vector<unsigned char>& mask,
                           int regionSize,
                           int paddedSize,
                           int maxLag,
                           const char* label) {
        std::vector<double> fftCovariance;
        if (!NoiseCovariance::computeMaskedCovarianceFFT(
                regionSize, paddedSize, maxLag, 0.25,
                residual, mask, fftCovariance)) {
            std::cerr << label << ": FFT covariance failed\n";
            return false;
        }
        const std::vector<double> direct = directCovariance(
            regionSize, maxLag, 0.25, residual, mask);
        return vectorsNear(fftCovariance, direct, 2.0e-11, label);
    }

    // ==========================================
    // Function: Exercise unmasked, random-mask, and production-like mask covariance tests
    // Method: Reuse a deterministic nontrivial field while progressively excluding pixels.
    // ==========================================
    bool testCovarianceAgainstDirect() {
        constexpr int regionSize = 12;
        constexpr int maxLag = 3;
        std::vector<double> residual(regionSize * regionSize, 0.0);
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                residual[static_cast<std::size_t>(y) * regionSize + x]
                    = std::sin(0.37 * x) + std::cos(0.23 * y) + 0.013 * x * y;
            }
        }

        std::vector<unsigned char> mask(regionSize * regionSize, 1U);
        if (!compareCovariance(residual, mask, regionSize, 2 * regionSize,
                               maxLag, "unmasked")) {
            return false;
        }

        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                if ((7 * x + 11 * y + 3) % 13 < 3) {
                    mask[static_cast<std::size_t>(y) * regionSize + x] = 0U;
                }
            }
        }
        if (!compareCovariance(residual, mask, regionSize, 2 * regionSize + 1,
                               maxLag, "random mask and rebuilt cache")) {
            return false;
        }

        std::fill(mask.begin(), mask.end(), 0U);
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize / 2; ++x) {
                const bool centralExclusion = x >= 3 && x < 6 && y >= 4 && y < 8;
                const bool dqHole = (x == 1 && (y == 2 || y == 9));
                if (!centralExclusion && !dqHole) {
                    mask[static_cast<std::size_t>(y) * regionSize + x] = 1U;
                }
            }
        }
        return compareCovariance(residual, mask, regionSize, 2 * regionSize,
                                 maxLag,
                                 "inner exclusion, DQ holes, and amplifier clipping");
    }

    // ==========================================
    // Function: Verify preservation of two-dimensional covariance anisotropy
    // Method: Use an x-alternating field and require horizontal/vertical lag estimates to differ.
    // ==========================================
    bool testAnisotropy() {
        constexpr int regionSize = 10;
        constexpr int maxLag = 2;
        std::vector<double> residual(regionSize * regionSize, 0.0);
        std::vector<unsigned char> mask(regionSize * regionSize, 1U);
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                residual[static_cast<std::size_t>(y) * regionSize + x]
                    = (x % 2 == 0 ? 1.0 : -1.0) + 0.02 * y;
            }
        }
        std::vector<double> covariance;
        if (!NoiseCovariance::computeMaskedCovarianceFFT(
                regionSize, 2 * regionSize, maxLag, 0.5,
                residual, mask, covariance)) {
            return false;
        }
        return std::abs(covariance[lagIndex(1, 0, maxLag)]
                        - covariance[lagIndex(0, 1, maxLag)]) > 0.5;
    }

    // ==========================================
    // Function: Verify covariance embedding, FFT normalization, and signed-power retention
    // Method: Compare a known short-lag covariance to its analytic shifted spectrum, then use a
    //         deliberately oscillatory covariance to require negative output modes.
    // ==========================================
    bool testCovarianceToPower() {
        constexpr int outputSize = 8;
        constexpr int maxLag = 1;
        std::vector<double> covariance(9, 0.0);
        covariance[lagIndex(0, 0, maxLag)] = 2.0;
        covariance[lagIndex(-1, 0, maxLag)] = 0.5;
        covariance[lagIndex(1, 0, maxLag)] = 0.5;
        covariance[lagIndex(0, -1, maxLag)] = 0.25;
        covariance[lagIndex(0, 1, maxLag)] = 0.25;

        std::vector<float> power;
        double maxImaginary = 0.0;
        double negativeFraction = 0.0;
        if (!NoiseCovariance::covarianceToNoisePower(
                outputSize, maxLag, covariance, power,
                maxImaginary, negativeFraction)) {
            return false;
        }
        if (maxImaginary > 1.0e-12 || negativeFraction != 0.0) {
            return false;
        }

        constexpr double pi = 3.14159265358979323846;
        for (int y = 0; y < outputSize; ++y) {
            const int kyIndex = (y + outputSize / 2) % outputSize;
            const double ky = 2.0 * pi * kyIndex / outputSize;
            for (int x = 0; x < outputSize; ++x) {
                const int kxIndex = (x + outputSize / 2) % outputSize;
                const double kx = 2.0 * pi * kxIndex / outputSize;
                const double expected = (2.0 + std::cos(kx) + 0.5 * std::cos(ky))
                                      / static_cast<double>(outputSize * outputSize);
                const double actual = power[static_cast<std::size_t>(y) * outputSize + x];
                if (std::abs(actual - expected) > 2.0e-8) {
                    std::cerr << "known covariance power mismatch\n";
                    return false;
                }
            }
        }

        std::fill(covariance.begin(), covariance.end(), 0.0);
        covariance[lagIndex(0, 0, maxLag)] = 0.1;
        covariance[lagIndex(-1, 0, maxLag)] = 1.0;
        covariance[lagIndex(1, 0, maxLag)] = 1.0;
        if (!NoiseCovariance::covarianceToNoisePower(
                outputSize, maxLag, covariance, power,
                maxImaginary, negativeFraction)) {
            return false;
        }
        return std::any_of(power.begin(), power.end(), [](float value) { return value < 0.0f; })
            && negativeFraction > 0.0 && maxImaginary <= 1.0e-12;
    }

    // ==========================================
    // Function: Guard the downstream stored-noise-power contract
    // Method: Copy a sentinel stamp containing negative modes and require exact preservation,
    //         proving that no FFT, magnitude square, absolute value, or clipping is applied.
    // ==========================================
    bool testStoredPowerCopy() {
        const std::vector<float> stored = {
            99.0f, 98.0f, 97.0f,
            -3.0f, 2.0f, -1.0f, 4.0f,
            96.0f
        };
        std::vector<float> copied;
        const std::vector<float> expected = {-3.0f, 2.0f, -1.0f, 4.0f};
        return NoiseCovariance::copyStoredNoisePower(stored, 3U, 2, copied)
            && copied == expected;
    }

    // ==========================================
    // Function: Verify the production 384-side padding contract
    // Method: Check the compile-time derivation, reject undersized padding, and compare the full
    //         192-side masked FFT covariance with direct pair summation at short lags.
    // ==========================================
    bool testProductionPadding() {
        static_assert(LensingConfig::noise_region_size == 192,
                      "production covariance region changed unexpectedly");
        static_assert(LensingConfig::noise_cov_padding_factor == 2.0,
                      "production covariance padding factor changed unexpectedly");
        static_assert(LensingConfig::noise_cov_fft_size == 384,
                      "production covariance FFT side must be 384");
        constexpr int regionSize = LensingConfig::noise_region_size;
        constexpr int maxLag = 2;
        std::vector<double> residual(regionSize * regionSize, 0.0);
        std::vector<unsigned char> mask(regionSize * regionSize, 1U);
        for (int y = 0; y < regionSize; ++y) {
            for (int x = 0; x < regionSize; ++x) {
                const std::size_t index = static_cast<std::size_t>(y) * regionSize + x;
                residual[index] = std::sin(0.031 * x) + std::cos(0.047 * y)
                                + 0.0003 * x * y;
                if ((5 * x + 7 * y) % 29 == 0) mask[index] = 0U;
            }
        }
        std::vector<double> rejected;
        if (NoiseCovariance::computeMaskedCovarianceFFT(
                regionSize, 2 * regionSize - 2, maxLag, 0.25,
                residual, mask, rejected)) {
            return false;
        }
        return compareCovariance(
            residual, mask, regionSize, LensingConfig::noise_cov_fft_size,
            maxLag, "production 384 padding");
    }

    // ==========================================
    // Function: Verify source-stamp amplifier-boundary rejection
    // Method: Exercise half-open intervals immediately left/right of the midpoint and one stamp
    //         that straddles it, while confirming single-amplifier mode remains unrestricted.
    // ==========================================
    bool testAmplifierBoundary() {
        constexpr int chipWidth = 2048;
        constexpr int stampSize = 64;
        constexpr int boundary = chipWidth / 2;
        return !NoiseCovariance::sourceStampCrossesAmplifier(
                   chipWidth, 2, boundary - stampSize, stampSize)
            && !NoiseCovariance::sourceStampCrossesAmplifier(
                   chipWidth, 2, boundary, stampSize)
            && NoiseCovariance::sourceStampCrossesAmplifier(
                   chipWidth, 2, boundary - stampSize / 2, stampSize)
            && !NoiseCovariance::sourceStampCrossesAmplifier(
                   chipWidth, 1, boundary - stampSize / 2, stampSize);
    }

}

// ==========================================
// Function: Run local-noise covariance and downstream-contract regression tests
// Method: Execute direct/FFT comparisons, mask geometry, anisotropy, analytic power, signed
//         mode, and verbatim stored-power checks; return nonzero on the first failed group.
// ==========================================
int main() {
    if (!testCovarianceAgainstDirect()) {
        std::cerr << "covariance comparison tests failed\n";
        return 1;
    }
    if (!testAnisotropy()) {
        std::cerr << "anisotropy test failed\n";
        return 1;
    }
    if (!testCovarianceToPower()) {
        std::cerr << "covariance-to-power tests failed\n";
        return 1;
    }
    if (!testStoredPowerCopy()) {
        std::cerr << "stored-power contract test failed\n";
        return 1;
    }
    if (!testProductionPadding()) {
        std::cerr << "production padding test failed\n";
        return 1;
    }
    if (!testAmplifierBoundary()) {
        std::cerr << "amplifier-boundary test failed\n";
        return 1;
    }
    std::cout << "Noise covariance tests passed\n";
    return 0;
}
