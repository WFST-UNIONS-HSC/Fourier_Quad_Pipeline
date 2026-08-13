#ifndef PSF_MODEL_STATE_HPP
#define PSF_MODEL_STATE_HPP

#include "LensingConfig.hpp"
#include "MPIFailure.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace PSFModel {
namespace Internal {

// ==========================================
// Structure: Store one chip's dynamically sized PSF candidates and chi grid
// Method: Reserve metadata initially, then allocate a live-stride row-major chi matrix.
// ==========================================
struct ChipPSFState {
    using StarRow = std::array<double, LensingConfig::npara>;

    std::vector<StarRow> stars;
    std::vector<float> chi_d;

    // ==========================================
    // Function: Allocate the full pairwise chi matrix for loaded candidates
    // Method: Guard size_t multiplication before zeroing one actual n-by-n matrix.
    // ==========================================
    void allocateChiD() {
        const std::size_t nstar = stars.size();
        if (nstar != 0
            && nstar > std::numeric_limits<std::size_t>::max() / nstar) {
            MPIFailure::abortWorld(
                "allocate PSF chi matrix", "candidate-count square overflow");
        }
        chi_d.assign(nstar * nstar, 0.0f);
    }

    // ==========================================
    // Function: Access one mutable pairwise chi value
    // Method: Use the live star count as the row-major matrix stride.
    // ==========================================
    float& getChiD(int star1, int star2) {
        const std::size_t nstar = stars.size();
        return chi_d[static_cast<std::size_t>(star1) * nstar
                     + static_cast<std::size_t>(star2)];
    }

    // ==========================================
    // Function: Access one immutable pairwise chi value
    // Method: Use the live star count as the row-major matrix stride.
    // ==========================================
    const float& getChiD(int star1, int star2) const {
        const std::size_t nstar = stars.size();
        return chi_d[static_cast<std::size_t>(star1) * nstar
                     + static_cast<std::size_t>(star2)];
    }
};

// ==========================================
// Structure: Own all per-chip PSF state for one exposure
// Method: Size only the live outer chip vector and derive star counts from dynamic storage.
// ==========================================
struct ExposurePSFState {
    std::vector<ChipPSFState> chips;

    // ==========================================
    // Function: Create empty dynamic state for the exposure's live chip count
    // Method: Avoid allocating any fixed per-chip candidate or chi matrix storage.
    // ==========================================
    explicit ExposurePSFState(int nchip)
        : chips(static_cast<std::size_t>(nchip)) {}

    // ==========================================
    // Function: Access one mutable star parameter
    // Method: Index the requested row in the chip's live candidate vector.
    // ==========================================
    double& getStarPara(int chip, int star, int para) {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Access one immutable star parameter
    // Method: Index the requested row in the chip's live candidate vector.
    // ==========================================
    const double& getStarPara(int chip, int star, int para) const {
        return chips[chip].stars[star][para];
    }

    // ==========================================
    // Function: Access one mutable exposure chi value
    // Method: Delegate to the chip's live-stride matrix accessor.
    // ==========================================
    float& getChiD(int chip, int star1, int star2) {
        return chips[chip].getChiD(star1, star2);
    }

    // ==========================================
    // Function: Access one immutable exposure chi value
    // Method: Delegate to the chip's live-stride matrix accessor.
    // ==========================================
    const float& getChiD(int chip, int star1, int star2) const {
        return chips[chip].getChiD(star1, star2);
    }

    // ==========================================
    // Function: Return one chip's live candidate count
    // Method: Derive the integer count from dynamic row storage.
    // ==========================================
    int getNStar(int chip) const {
        return static_cast<int>(chips[chip].stars.size());
    }
};

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_MODEL_STATE_HPP
