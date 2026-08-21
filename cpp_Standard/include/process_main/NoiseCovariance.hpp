#ifndef NOISE_COVARIANCE_HPP
#define NOISE_COVARIANCE_HPP

#include <cstddef>
#include <vector>

namespace NoiseCovariance {

    bool computeMaskedCovarianceFFT(
        int regionSize, int paddedSize, int maxLag, double minPairFraction,
        const std::vector<double>& residual,
        const std::vector<unsigned char>& mask,
        std::vector<double>& covariance,
        std::vector<double>* pairCounts = nullptr);

    // ==========================================
    // Function: Detect a source stamp that straddles a two-amplifier boundary
    // Method: Compare the half-open stamp x interval with the chip midpoint only when split mode 2
    //         is active; single-amplifier mode never rejects on this geometry.
    // ==========================================
    bool sourceStampCrossesAmplifier(
        int chipWidth, int ccdSplit, int stampStartX, int stampSize);

    bool covarianceToNoisePower(
        int outputSize, int maxLag,
        const std::vector<double>& covariance,
        std::vector<float>& noisePower,
        double& maxImaginary,
        double& negativeFraction);

    bool copyStoredNoisePower(
        const std::vector<float>& storedNoisePower,
        std::size_t offset,
        int stampSize,
        std::vector<float>& noisePower);

}

#endif // NOISE_COVARIANCE_HPP
