#ifndef EXPOSURE_INFO_HPP
#define EXPOSURE_INFO_HPP

#include <cstddef>
#include <vector>
#include <string>

namespace ExposureInfo {
    struct State {
        std::vector<float> parameters;

        // ==========================================
        // Function: Reset Stage-8 aggregate storage
        // Method: Allocate exactly the requested number of zeroed values.
        // ==========================================
        void reset(std::size_t count) {
            parameters.assign(count, 0.0f);
        }
    };

    extern State state;

    // ==========================================
    // Function: Return the live Stage-8 parameter count
    // Method: Allocate six aggregate values per runtime exposure and none for invalid counts.
    // ==========================================
    inline std::size_t parameterCount(int exposure_count) {
        return exposure_count > 0
            ? static_cast<std::size_t>(exposure_count) * 6
            : 0;
    }

    void getExpoInfo(const std::vector<std::string>& imageFiles, int nchip, const std::string& dirOutput, float para[6]);
    void procInfo(int iexpo);
}

#endif // EXPOSURE_INFO_HPP
