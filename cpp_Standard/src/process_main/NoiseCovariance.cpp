#include "NoiseCovariance.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

#include <fftw3.h>

namespace NoiseCovariance {

    namespace {

        // ==========================================
        // Class: Reusable zero-padded autocorrelation workspace
        // Method: Bind fixed-size FFTW plans to two complex buffers and reuse them for every
        //         residual/mask pair evaluated by one execution thread.
        // ==========================================
        class AutocorrelationWorkspace {
        public:
            // ==========================================
            // Function: Construct fixed-size forward/backward autocorrelation plans
            // Method: Allocate full linear-correlation padding of side 2*N-1 and bind in-place
            //         FFTW transforms for the residual and binary mask buffers.
            // ==========================================
            explicit AutocorrelationWorkspace(int regionSize)
                : regionSize_(regionSize),
                  paddedSize_(2 * regionSize - 1),
                  residualBuffer_(elementCount()),
                  maskBuffer_(elementCount()) {
                residualForward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, residualData(), residualData(),
                    FFTW_FORWARD, FFTW_ESTIMATE);
                residualBackward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, residualData(), residualData(),
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                maskForward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, maskData(), maskData(),
                    FFTW_FORWARD, FFTW_ESTIMATE);
                maskBackward_ = fftw_plan_dft_2d(
                    paddedSize_, paddedSize_, maskData(), maskData(),
                    FFTW_BACKWARD, FFTW_ESTIMATE);
            }

            // ==========================================
            // Function: Destroy reusable autocorrelation plans
            // Method: Release every successfully constructed FFTW plan while vector storage is
            //         still alive.
            // ==========================================
            ~AutocorrelationWorkspace() {
                if (residualForward_ != nullptr) fftw_destroy_plan(residualForward_);
                if (residualBackward_ != nullptr) fftw_destroy_plan(residualBackward_);
                if (maskForward_ != nullptr) fftw_destroy_plan(maskForward_);
                if (maskBackward_ != nullptr) fftw_destroy_plan(maskBackward_);
            }

            AutocorrelationWorkspace(const AutocorrelationWorkspace&) = delete;
            AutocorrelationWorkspace& operator=(const AutocorrelationWorkspace&) = delete;

            // ==========================================
            // Function: Report whether all autocorrelation plans were created
            // Method: Require the four bound FFTW plan handles used by execute().
            // ==========================================
            bool valid() const {
                return residualForward_ != nullptr && residualBackward_ != nullptr
                    && maskForward_ != nullptr && maskBackward_ != nullptr;
            }

            // ==========================================
            // Function: Compute linear residual and mask autocorrelations
            // Method: Zero-pad at the origin, transform both fields, replace each spectrum by
            //         its squared magnitude, and inverse transform in place.
            // ==========================================
            void execute(const std::vector<double>& residual,
                         const std::vector<unsigned char>& mask) {
                std::fill(residualBuffer_.begin(), residualBuffer_.end(),
                          std::complex<double>(0.0, 0.0));
                std::fill(maskBuffer_.begin(), maskBuffer_.end(),
                          std::complex<double>(0.0, 0.0));

                for (int y = 0; y < regionSize_; ++y) {
                    for (int x = 0; x < regionSize_; ++x) {
                        const std::size_t inputIndex = static_cast<std::size_t>(y)
                                                     * static_cast<std::size_t>(regionSize_)
                                                     + static_cast<std::size_t>(x);
                        const std::size_t paddedIndex = static_cast<std::size_t>(y)
                                                      * static_cast<std::size_t>(paddedSize_)
                                                      + static_cast<std::size_t>(x);
                        if (mask[inputIndex] != 0U) {
                            residualBuffer_[paddedIndex] = residual[inputIndex];
                            maskBuffer_[paddedIndex] = 1.0;
                        }
                    }
                }

                fftw_execute(residualForward_);
                fftw_execute(maskForward_);
                for (std::size_t i = 0; i < elementCount(); ++i) {
                    residualBuffer_[i] = std::complex<double>(
                        std::norm(residualBuffer_[i]), 0.0);
                    maskBuffer_[i] = std::complex<double>(
                        std::norm(maskBuffer_[i]), 0.0);
                }
                fftw_execute(residualBackward_);
                fftw_execute(maskBackward_);
            }

            // ==========================================
            // Function: Read one normalized autocorrelation numerator
            // Method: Map a signed linear lag to the periodic zero-padded FFT index and remove
            //         FFTW's unnormalized inverse-transform scale.
            // ==========================================
            double residualCorrelation(int dx, int dy) const {
                return residualBuffer_[lagIndex(dx, dy)].real()
                     / static_cast<double>(elementCount());
            }

            // ==========================================
            // Function: Read one normalized valid-pair count
            // Method: Apply the same lag mapping and inverse scale as the residual correlation.
            // ==========================================
            double maskCorrelation(int dx, int dy) const {
                return maskBuffer_[lagIndex(dx, dy)].real()
                     / static_cast<double>(elementCount());
            }

            // ==========================================
            // Function: Return the configured source-region side length
            // Method: Expose the immutable plan key for thread-local cache reuse.
            // ==========================================
            int regionSize() const {
                return regionSize_;
            }

        private:
            // ==========================================
            // Function: Count complex elements in one padded FFT buffer
            // Method: Square the fixed padded side in size_t arithmetic.
            // ==========================================
            std::size_t elementCount() const {
                return static_cast<std::size_t>(paddedSize_)
                     * static_cast<std::size_t>(paddedSize_);
            }

            // ==========================================
            // Function: Convert a signed lag to a padded row-major index
            // Method: Wrap negative x/y lags to the upper periodic indices of the full
            //         2*N-1 linear-correlation grid.
            // ==========================================
            std::size_t lagIndex(int dx, int dy) const {
                const int x = dx >= 0 ? dx : paddedSize_ + dx;
                const int y = dy >= 0 ? dy : paddedSize_ + dy;
                return static_cast<std::size_t>(y) * static_cast<std::size_t>(paddedSize_)
                     + static_cast<std::size_t>(x);
            }

            // ==========================================
            // Function: Access residual storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* residualData() {
                return reinterpret_cast<fftw_complex*>(residualBuffer_.data());
            }

            // ==========================================
            // Function: Access mask storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* maskData() {
                return reinterpret_cast<fftw_complex*>(maskBuffer_.data());
            }

            int regionSize_ = 0;
            int paddedSize_ = 0;
            std::vector<std::complex<double>> residualBuffer_;
            std::vector<std::complex<double>> maskBuffer_;
            fftw_plan residualForward_ = nullptr;
            fftw_plan residualBackward_ = nullptr;
            fftw_plan maskForward_ = nullptr;
            fftw_plan maskBackward_ = nullptr;
        };

        // ==========================================
        // Class: Reusable covariance-to-power workspace
        // Method: Bind one fixed-size forward FFTW plan to the output covariance grid.
        // ==========================================
        class PowerWorkspace {
        public:
            // ==========================================
            // Function: Construct a fixed-size covariance transform
            // Method: Allocate one square complex grid and bind an in-place forward FFTW plan.
            // ==========================================
            explicit PowerWorkspace(int outputSize)
                : outputSize_(outputSize),
                  buffer_(static_cast<std::size_t>(outputSize)
                          * static_cast<std::size_t>(outputSize)) {
                forward_ = fftw_plan_dft_2d(
                    outputSize_, outputSize_, data(), data(), FFTW_FORWARD, FFTW_ESTIMATE);
            }

            // ==========================================
            // Function: Destroy the reusable covariance transform
            // Method: Release the FFTW plan before its bound vector storage is destroyed.
            // ==========================================
            ~PowerWorkspace() {
                if (forward_ != nullptr) fftw_destroy_plan(forward_);
            }

            PowerWorkspace(const PowerWorkspace&) = delete;
            PowerWorkspace& operator=(const PowerWorkspace&) = delete;

            // ==========================================
            // Function: Report whether the covariance transform plan exists
            // Method: Validate the one forward FFTW plan required by execute().
            // ==========================================
            bool valid() const {
                return forward_ != nullptr;
            }

            // ==========================================
            // Function: Transform an unshifted covariance grid
            // Method: Copy the caller grid into bound storage and execute the reusable forward
            //         transform without altering signed Fourier modes.
            // ==========================================
            void execute(const std::vector<std::complex<double>>& covarianceGrid) {
                std::copy(covarianceGrid.begin(), covarianceGrid.end(), buffer_.begin());
                fftw_execute(forward_);
            }

            // ==========================================
            // Function: Access one transformed Fourier coefficient
            // Method: Read the immutable result buffer by row-major index.
            // ==========================================
            const std::complex<double>& value(std::size_t index) const {
                return buffer_[index];
            }

            // ==========================================
            // Function: Return the configured Fourier-grid side length
            // Method: Expose the immutable plan key for thread-local cache reuse.
            // ==========================================
            int outputSize() const {
                return outputSize_;
            }

        private:
            // ==========================================
            // Function: Access power-grid storage through FFTW's complex ABI
            // Method: Reinterpret the bound std::complex buffer used throughout this codebase.
            // ==========================================
            fftw_complex* data() {
                return reinterpret_cast<fftw_complex*>(buffer_.data());
            }

            int outputSize_ = 0;
            std::vector<std::complex<double>> buffer_;
            fftw_plan forward_ = nullptr;
        };

        // ==========================================
        // Function: Reuse the current thread's autocorrelation workspace
        // Method: Rebuild the cached workspace only when the configured region size changes.
        // ==========================================
        AutocorrelationWorkspace& autocorrelationWorkspace(int regionSize) {
            thread_local std::unique_ptr<AutocorrelationWorkspace> workspace;
            if (!workspace || workspace->regionSize() != regionSize) {
                workspace = std::make_unique<AutocorrelationWorkspace>(regionSize);
            }
            return *workspace;
        }

        // ==========================================
        // Function: Reuse the current thread's covariance-to-power workspace
        // Method: Rebuild the cached workspace only when the output stamp size changes.
        // ==========================================
        PowerWorkspace& powerWorkspace(int outputSize) {
            thread_local std::unique_ptr<PowerWorkspace> workspace;
            if (!workspace || workspace->outputSize() != outputSize) {
                workspace = std::make_unique<PowerWorkspace>(outputSize);
            }
            return *workspace;
        }

        // ==========================================
        // Function: Map a signed lag to the compact covariance-vector index
        // Method: Offset both axes by maxLag in a (2*maxLag+1)^2 row-major grid.
        // ==========================================
        std::size_t compactLagIndex(int dx, int dy, int maxLag) {
            const int side = 2 * maxLag + 1;
            return static_cast<std::size_t>(dy + maxLag) * static_cast<std::size_t>(side)
                 + static_cast<std::size_t>(dx + maxLag);
        }

    }

    // ==========================================
    // Function: Estimate a masked two-dimensional covariance with FFT acceleration
    // Method: Compute full zero-padded linear autocorrelations for M*R and M, divide by valid
    //         pair counts, retain qualified short lags, and enforce C(d)=C(-d) symmetry.
    // ==========================================
    bool computeMaskedCovarianceFFT(
        int regionSize, int maxLag, double minPairFraction,
        const std::vector<double>& residual,
        const std::vector<unsigned char>& mask,
        std::vector<double>& covariance,
        std::vector<double>* pairCounts) {
        if (regionSize <= 0 || maxLag < 0 || maxLag >= regionSize
            || !std::isfinite(minPairFraction)
            || minPairFraction <= 0.0 || minPairFraction > 1.0) {
            return false;
        }

        const std::size_t regionElements = static_cast<std::size_t>(regionSize)
                                         * static_cast<std::size_t>(regionSize);
        if (residual.size() != regionElements || mask.size() != regionElements) {
            return false;
        }
        for (std::size_t i = 0; i < regionElements; ++i) {
            if (mask[i] != 0U && !std::isfinite(residual[i])) {
                return false;
            }
        }

        AutocorrelationWorkspace& workspace = autocorrelationWorkspace(regionSize);
        if (!workspace.valid()) {
            return false;
        }
        workspace.execute(residual, mask);

        const int lagSide = 2 * maxLag + 1;
        const std::size_t lagElements = static_cast<std::size_t>(lagSide)
                                      * static_cast<std::size_t>(lagSide);
        covariance.assign(lagElements, 0.0);
        std::vector<double> localPairCounts(lagElements, 0.0);

        const double zeroLagPairs = workspace.maskCorrelation(0, 0);
        if (!std::isfinite(zeroLagPairs) || zeroLagPairs < 1.0) {
            return false;
        }
        const double minimumPairs = minPairFraction * zeroLagPairs;

        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                const std::size_t index = compactLagIndex(dx, dy, maxLag);
                const double pairs = std::max(0.0, workspace.maskCorrelation(dx, dy));
                localPairCounts[index] = pairs;
                if (pairs + 1.0e-8 < minimumPairs) {
                    continue;
                }
                const double numerator = workspace.residualCorrelation(dx, dy);
                const double value = numerator / pairs;
                if (!std::isfinite(value)) {
                    return false;
                }
                covariance[index] = value;
            }
        }

        for (int dy = 0; dy <= maxLag; ++dy) {
            const int dxStart = dy == 0 ? 0 : -maxLag;
            for (int dx = dxStart; dx <= maxLag; ++dx) {
                const std::size_t forwardIndex = compactLagIndex(dx, dy, maxLag);
                const std::size_t reverseIndex = compactLagIndex(-dx, -dy, maxLag);
                const bool forwardValid = localPairCounts[forwardIndex] + 1.0e-8 >= minimumPairs;
                const bool reverseValid = localPairCounts[reverseIndex] + 1.0e-8 >= minimumPairs;
                if (forwardValid && reverseValid) {
                    const double symmetricValue = 0.5 * (
                        covariance[forwardIndex] + covariance[reverseIndex]);
                    covariance[forwardIndex] = symmetricValue;
                    covariance[reverseIndex] = symmetricValue;
                } else {
                    covariance[forwardIndex] = 0.0;
                    covariance[reverseIndex] = 0.0;
                }
            }
        }

        if (pairCounts != nullptr) {
            *pairCounts = std::move(localPairCounts);
        }
        return true;
    }

    // ==========================================
    // Function: Convert a short-lag covariance estimate into stored noise power
    // Method: Embed signed lags at periodic indices, apply a normalized forward FFT matching
    //         getPower(), fftshift the real result, and report imaginary/negative diagnostics.
    // ==========================================
    bool covarianceToNoisePower(
        int outputSize, int maxLag,
        const std::vector<double>& covariance,
        std::vector<float>& noisePower,
        double& maxImaginary,
        double& negativeFraction) {
        if (outputSize <= 0 || maxLag < 0 || maxLag >= outputSize / 2) {
            return false;
        }
        const int lagSide = 2 * maxLag + 1;
        const std::size_t expectedCovariance = static_cast<std::size_t>(lagSide)
                                               * static_cast<std::size_t>(lagSide);
        if (covariance.size() != expectedCovariance) {
            return false;
        }

        const std::size_t outputElements = static_cast<std::size_t>(outputSize)
                                         * static_cast<std::size_t>(outputSize);
        std::vector<std::complex<double>> covarianceGrid(
            outputElements, std::complex<double>(0.0, 0.0));
        for (int dy = -maxLag; dy <= maxLag; ++dy) {
            const int y = dy >= 0 ? dy : outputSize + dy;
            for (int dx = -maxLag; dx <= maxLag; ++dx) {
                const int x = dx >= 0 ? dx : outputSize + dx;
                const double value = covariance[compactLagIndex(dx, dy, maxLag)];
                if (!std::isfinite(value)) {
                    return false;
                }
                covarianceGrid[static_cast<std::size_t>(y)
                               * static_cast<std::size_t>(outputSize)
                               + static_cast<std::size_t>(x)] = value;
            }
        }

        PowerWorkspace& workspace = powerWorkspace(outputSize);
        if (!workspace.valid()) {
            return false;
        }
        workspace.execute(covarianceGrid);

        noisePower.assign(outputElements, 0.0f);
        maxImaginary = 0.0;
        negativeFraction = 0.0;
        double negativeMagnitude = 0.0;
        double totalMagnitude = 0.0;
        const double normalization = 1.0 / static_cast<double>(outputElements);
        const int half = outputSize / 2;

        for (int y = 0; y < outputSize; ++y) {
            const int shiftedY = (y + half) % outputSize;
            for (int x = 0; x < outputSize; ++x) {
                const int shiftedX = (x + half) % outputSize;
                const std::size_t inputIndex = static_cast<std::size_t>(y)
                                             * static_cast<std::size_t>(outputSize)
                                             + static_cast<std::size_t>(x);
                const std::size_t outputIndex = static_cast<std::size_t>(shiftedY)
                                              * static_cast<std::size_t>(outputSize)
                                              + static_cast<std::size_t>(shiftedX);
                const std::complex<double> value = workspace.value(inputIndex) * normalization;
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                    return false;
                }
                maxImaginary = std::max(maxImaginary, std::abs(value.imag()));
                totalMagnitude += std::abs(value.real());
                if (value.real() < 0.0) {
                    negativeMagnitude += -value.real();
                }
                noisePower[outputIndex] = static_cast<float>(value.real());
            }
        }

        if (totalMagnitude > 0.0) {
            negativeFraction = negativeMagnitude / totalMagnitude;
        }
        return std::isfinite(negativeFraction);
    }

    // ==========================================
    // Function: Load one stored signed noise-power stamp without re-transforming it
    // Method: Validate the requested contiguous ns-by-ns slice and copy its values verbatim,
    //         preserving negative Fourier modes for linear subtraction.
    // ==========================================
    bool copyStoredNoisePower(
        const std::vector<float>& storedNoisePower,
        std::size_t offset,
        int stampSize,
        std::vector<float>& noisePower) {
        if (stampSize <= 0) {
            return false;
        }
        const std::size_t stampElements = static_cast<std::size_t>(stampSize)
                                        * static_cast<std::size_t>(stampSize);
        if (offset > storedNoisePower.size()
            || stampElements > storedNoisePower.size() - offset) {
            return false;
        }
        noisePower.assign(storedNoisePower.begin() + static_cast<std::ptrdiff_t>(offset),
                          storedNoisePower.begin()
                              + static_cast<std::ptrdiff_t>(offset + stampElements));
        return true;
    }

}
