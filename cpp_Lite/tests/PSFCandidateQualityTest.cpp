#include "process_main/PSFCandidateQuality.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using namespace PSFModel::Internal;

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF candidate quality test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void testF77CompatibilityPolicy() {
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    std::vector<float> power(25, 0.0f);
    power[12] = 20.0f;
    const int neighbors[] = {6, 7, 8, 11, 13, 16, 17, 18};
    for (int index : neighbors) power[index] = -1.0f;
    require(assessF77CandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::Accepted,
            "negative core neighbours must not trigger a modern admission gate");

    power[0] = -20.0f;
    require(assessF77CandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::Accepted
                && sum_power < 0.0 && chi_window_sum < 0.0,
            "finite negative sums must remain normalizable in F77 mode");

    power[0] = -12.0f;
    require(assessF77CandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonPositiveSum
                && sum_power == 0.0,
            "a zero full sum must fail the normalization guard");

    power[0] = std::numeric_limits<float>::quiet_NaN();
    require(assessF77CandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonFinitePower,
            "non-finite power must be rejected");
}

void testStructuralAndDiagnosticValidity() {
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    require(assessF77CandidatePower(
                5, 5, std::vector<float>(24), sum_power, chi_window_sum)
                == CandidatePowerStatus::InvalidShape,
            "incomplete power image must be rejected");
    require(assessF77CandidatePower(
                2, 2, std::vector<float>(4), sum_power, chi_window_sum)
                == CandidatePowerStatus::InvalidShape,
            "too-small image must be rejected");
    require(candidateDiagnosticsAreFinite(0.0, 0.1, -0.1),
            "finite zero size must not be rejected by a modern positivity gate");
    require(!candidateDiagnosticsAreFinite(
                1.0, std::numeric_limits<double>::quiet_NaN(), 0.0),
            "NaN ellipticity must be rejected");
    require(!candidateDiagnosticsAreFinite(
                std::numeric_limits<double>::infinity(), 0.0, 0.0),
            "infinite size must be rejected");
}

}  // namespace

int main() {
    testF77CompatibilityPolicy();
    testStructuralAndDiagnosticValidity();
    std::cout << "PSF candidate quality tests passed\n";
    return EXIT_SUCCESS;
}
