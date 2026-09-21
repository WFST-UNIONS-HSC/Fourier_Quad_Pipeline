#include "process_main/PSFModelState.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSFModelState test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void testChipSize(int star_count) {
    PSFModel::Internal::ExposurePSFState state(1);
    auto& chip = state.chips[0];
    chip.stars.reserve(LensingConfig::nstar_max);
    for (int star = 0; star < star_count; ++star) {
        PSFModel::Internal::ChipPSFState::StarRow row{};
        row[0] = star;
        chip.stars.push_back(row);
    }
    chip.selection.resize(static_cast<std::size_t>(star_count));
    for (int star = 0; star < star_count; ++star) {
        chip.selection[star].chi_window.push_back(static_cast<float>(star));
    }
    require(state.getNStar(0) == star_count,
            "candidate count must follow dynamic storage");
    require(chip.selection.size() == static_cast<std::size_t>(star_count),
            "selection metadata must align with candidate rows");
    if (star_count > 0) {
        require(chip.selection.back().chi_window.size() == 1
                    && !chip.selection.back().selected_group
                    && !chip.selection.back().selected_fit
                    && chip.selection.back().leverage == 0.0,
                "candidate state must default to the frozen selection policy");
    }
}

void testDynamicChipSizes() {
    const int counts[] = {0, 10, 1999, 2000, 2001, 2301};
    for (int count : counts) testChipSize(count);
}

void testInitialFitCacheClear() {
    PSFModel::Internal::ChipPSFFitState fit;
    fit.valid = true;
    fit.initial_star_count = 3;
    fit.star_indices = {1, 2, 3};
    fit.coefficients = {10.0, 20.0};
    fit.leverage = {0.1, 0.2, 0.3};
    fit.clear();
    require(!fit.valid && fit.initial_star_count == 0
                && fit.star_indices.empty() && fit.coefficients.empty()
                && fit.leverage.empty(),
            "clearing a failed initial fit must remove only cached fit products");
}

}  // namespace

int main() {
    testDynamicChipSizes();
    testInitialFitCacheClear();
    std::cout << "PSFModelState tests passed\n";
    return EXIT_SUCCESS;
}
