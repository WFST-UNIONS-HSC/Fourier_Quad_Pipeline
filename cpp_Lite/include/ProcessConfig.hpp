#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include <string>
#include <vector>

namespace ProcessConfig {

struct DatasetSpec {
    std::string target;
    std::string prefix;
};

// Workflow defaults. Command-line --run-init/--run-main values override these.
inline constexpr bool RUN_PROCESS_INIT = false;
inline constexpr bool RUN_PROCESS_MAIN = true;

// Initializer and exposure-list defaults. Edit this file for one site's usual dataset.
inline constexpr const char* SCIENCE_ROOT = "/lustre/home/acct-phyzj/share/DES/g";
inline constexpr const char* DQ_ROOT = "/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask";
inline constexpr const char* OUTPUT_ROOT = "/lustre/home/acct-phyzj/share/DES/g_band_v1";
inline const std::vector<DatasetSpec> DATASETS = {{"g2019", "c4d_19"}};
inline const std::vector<std::string> CONTAINS = {"v1"};
inline constexpr const char* EXISTING = "fail";
inline constexpr int F77_MAX_PATH = 149;
inline constexpr const char* EXPO_LIST = "";

// ==========================================
// Configuration: Runtime workflow options
// Method: Seed every optional command-line value from the defaults above, then
//         let the unified parser override only explicitly supplied options.
// ==========================================
struct RuntimeOptions {
    bool run_process_init = RUN_PROCESS_INIT;
    bool run_process_main = RUN_PROCESS_MAIN;
    std::string science_root = SCIENCE_ROOT;
    std::string dq_root = DQ_ROOT;
    std::string output_root = OUTPUT_ROOT;
    std::vector<DatasetSpec> datasets = DATASETS;
    std::vector<std::string> contains = CONTAINS;
    std::string existing = EXISTING;
    int f77_max_path = F77_MAX_PATH;
    std::string expo_list = EXPO_LIST;
    bool external_expo_list_supplied = false;
    bool help_requested = false;
};

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
