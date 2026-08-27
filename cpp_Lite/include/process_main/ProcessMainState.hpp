#ifndef PROCESS_MAIN_STATE_HPP
#define PROCESS_MAIN_STATE_HPP

#include <string>
#include <vector>

namespace ProcessMain {

struct State {
    std::vector<std::string> exposure_files;
    std::string source_catalog_directory;

    // ==========================================
    // Function: Return the live exposure count
    // Method: Derive the count from the owning vector to avoid duplicate state.
    // ==========================================
    int exposureCount() const noexcept {
        return static_cast<int>(exposure_files.size());
    }

    // ==========================================
    // Function: Clear per-run process_main state
    // Method: Release exposure paths and the effective catalog directory
    //         together before a new dataset begins.
    // ==========================================
    void clear() {
        exposure_files.clear();
        source_catalog_directory.clear();
    }
};

extern State state;

}  // namespace ProcessMain

#endif  // PROCESS_MAIN_STATE_HPP
