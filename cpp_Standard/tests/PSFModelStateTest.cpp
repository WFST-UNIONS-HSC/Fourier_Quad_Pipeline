#include "PSFModelState.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

// ==========================================
// Function: Stop the PSF-state test when one invariant fails
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFModelState test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Verify one actual candidate count and chi-matrix dimension
// Method: Reserve the legacy hint, append all rows, and exercise the last live-stride element.
// ==========================================
void testChipSize(int star_count) {
    PSFModel::Internal::ExposurePSFState state(1);
    auto& chip = state.chips[0];
    chip.stars.reserve(LensingConfig::nstar_max);
    for (int star = 0; star < star_count; ++star) {
        PSFModel::Internal::ChipPSFState::StarRow row{};
        row[0] = star;
        chip.stars.push_back(row);
    }
    chip.allocateChiD();

    const std::size_t expected =
        static_cast<std::size_t>(star_count) * star_count;
    require(state.getNStar(0) == star_count,
            "candidate count must follow dynamic storage");
    require(chip.chi_d.size() == expected,
            "chi matrix must use actual nstar squared");
    if (star_count > 1) {
        state.getChiD(0, star_count - 1, star_count - 2) = 7.5f;
        require(state.getChiD(0, star_count - 1, star_count - 2) == 7.5f,
                "last live-stride chi element must be writable");
    }
}

// ==========================================
// Function: Run focused Stage-5 dynamic-state regression cases
// Method: Exercise zero, ordinary, boundary, and above-reservation star counts.
// ==========================================
void testDynamicChipSizes() {
    const int counts[] = {0, 10, 1999, 2000, 2001, 2301};
    for (int count : counts) {
        testChipSize(count);
    }
}

}  // namespace

// ==========================================
// Function: Run the PSF-state test suite
// Method: Execute all live-size cases and report one success line.
// ==========================================
int main() {
    testDynamicChipSizes();
    std::cout << "PSFModelState tests passed\n";
    return EXIT_SUCCESS;
}
