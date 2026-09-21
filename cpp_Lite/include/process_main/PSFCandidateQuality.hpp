#ifndef PSF_CANDIDATE_QUALITY_HPP
#define PSF_CANDIDATE_QUALITY_HPP

#include "process_main/PSFStarSelection.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace PSFModel {
namespace Internal {

enum class CandidatePowerStatus {
    Accepted,
    InvalidShape,
    NonFinitePower,
    NonPositiveSum
};

// ==========================================
// Function: Apply only hard numerical checks needed by F77 PSF selection
// Method: Require a finite square stamp and finite nonzero full sum while
//         deliberately omitting modern core-median and chi-window-sign gates.
// ==========================================
inline CandidatePowerStatus assessF77CandidatePower(
    int nx, int ny, const std::vector<float>& power,
    double& sum_power, double& chi_window_sum) {
    sum_power = 0.0;
    chi_window_sum = 0.0;
    if (nx < 3 || ny < 3 || nx != ny
        || power.size() != static_cast<std::size_t>(nx) * ny) {
        return CandidatePowerStatus::InvalidShape;
    }
    for (float value : power) {
        if (!std::isfinite(value)) return CandidatePowerStatus::NonFinitePower;
        sum_power += static_cast<double>(value);
    }
    if (!std::isfinite(sum_power) || sum_power == 0.0) {
        return CandidatePowerStatus::NonPositiveSum;
    }
    const PSFChiWindow chi_window = getPSFChiWindow(nx);
    if (chi_window.pixelCount() <= 0) {
        return CandidatePowerStatus::InvalidShape;
    }
    for (int row = chi_window.first; row <= chi_window.last; ++row) {
        for (int column = chi_window.first; column <= chi_window.last; ++column) {
            chi_window_sum += static_cast<double>(power[row * nx + column]);
        }
    }
    return std::isfinite(chi_window_sum)
        ? CandidatePowerStatus::Accepted
        : CandidatePowerStatus::NonFinitePower;
}

// ==========================================
// Function: Validate diagnostics derived from an accepted candidate spectrum
// Method: Reject only non-finite results so unsafe values cannot reach grouping.
// ==========================================
inline bool candidateDiagnosticsAreFinite(
    double size, double e1, double e2) {
    return std::isfinite(size) && std::isfinite(e1)
        && std::isfinite(e2);
}

}  // namespace Internal
}  // namespace PSFModel

#endif  // PSF_CANDIDATE_QUALITY_HPP
