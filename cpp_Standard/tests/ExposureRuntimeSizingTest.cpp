#include "process_main/ExposureInfo.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace ExposureInfo {
State state;
}

namespace {

// ==========================================
// Function: Stop the Stage-8 sizing test on a failed requirement
// Method: Print one focused diagnostic and terminate with failure status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "Exposure runtime-sizing test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Verify one runtime exposure count
// Method: Use the production count helper for allocation and require six slots per exposure.
// ==========================================
void testCount(int exposure_count) {
    const std::size_t expected = static_cast<std::size_t>(exposure_count) * 6;
    const std::size_t count = ExposureInfo::parameterCount(exposure_count);
    ExposureInfo::state.parameters.assign(count, 0.0f);
    require(count == expected && ExposureInfo::state.parameters.size() == expected,
            "unexpected Stage-8 vector or MPI element count");
}

}  // namespace

// ==========================================
// Function: Run Stage-8 runtime sizing regressions
// Method: Cover representative small and multi-exposure production counts.
// ==========================================
int main() {
    require(ExposureInfo::parameterCount(0) == 0,
            "zero exposures must allocate no Stage-8 storage");
    const int exposure_counts[] = {1, 2, 17, 100};
    for (int exposure_count : exposure_counts) {
        testCount(exposure_count);
    }
    std::cout << "Exposure runtime-sizing tests passed\n";
    return EXIT_SUCCESS;
}
