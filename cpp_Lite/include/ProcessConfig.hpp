#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include "process_main/LensingConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ProcessConfig {

struct DatasetSpec {
    std::string target;
    std::string prefix;
};

// Workflow defaults. Command-line phase switches override these values.
inline constexpr bool RUN_PROCESS_EXTCAT = false;
inline constexpr bool RUN_PROCESS_INIT = false;
inline constexpr bool RUN_PROCESS_MAIN = true;

// ==========================================
// Configuration: External source-catalog repartitioning defaults
// Method: Derive the tile output from the main pipeline SOURCE_CAT while keeping raw-input,
//         parsing, and optional projection controls together for the first workflow phase.
// ==========================================
inline constexpr const char* EXTCAT_INPUT_DIRECTORY = "";
inline const std::string& EXTCAT_OUTPUT_DIRECTORY = LensingConfig::SOURCE_CAT;
inline const std::vector<std::string> EXTCAT_FILENAME_TOKENS = {};
inline constexpr bool EXTCAT_RECURSIVE = true;
inline constexpr const char* EXTCAT_DELIMITER = "auto";
inline constexpr const char* EXTCAT_HEADER_MODE = "auto";
inline constexpr const char* EXTCAT_MALFORMED_POLICY = "fail";
inline constexpr const char* EXTCAT_EXISTING_POLICY = "fail";
inline constexpr std::uint64_t EXTCAT_CHUNK_MIB = 64;
inline constexpr bool EXTCAT_USE_EXPLICIT_COLUMNS = false;
inline const std::vector<std::size_t> EXTCAT_INPUT_COLUMNS_ONE_BASED = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
};
inline constexpr bool EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS = false;
inline constexpr std::size_t EXTCAT_RA_COLUMN_ONE_BASED = 5;
inline constexpr std::size_t EXTCAT_DEC_COLUMN_ONE_BASED = 6;

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
    bool run_process_extcat = RUN_PROCESS_EXTCAT;
    bool run_process_init = RUN_PROCESS_INIT;
    bool run_process_main = RUN_PROCESS_MAIN;
    std::string extcat_input_directory = EXTCAT_INPUT_DIRECTORY;
    std::string extcat_output_directory = EXTCAT_OUTPUT_DIRECTORY;
    std::vector<std::string> extcat_filename_tokens = EXTCAT_FILENAME_TOKENS;
    bool extcat_recursive = EXTCAT_RECURSIVE;
    std::string extcat_delimiter = EXTCAT_DELIMITER;
    std::string extcat_header_mode = EXTCAT_HEADER_MODE;
    std::string extcat_malformed_policy = EXTCAT_MALFORMED_POLICY;
    std::string extcat_existing_policy = EXTCAT_EXISTING_POLICY;
    std::uint64_t extcat_chunk_mib = EXTCAT_CHUNK_MIB;
    bool extcat_use_explicit_columns = EXTCAT_USE_EXPLICIT_COLUMNS;
    std::vector<std::size_t> extcat_input_columns_one_based =
        EXTCAT_INPUT_COLUMNS_ONE_BASED;
    bool extcat_use_explicit_coordinate_columns =
        EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS;
    std::size_t extcat_ra_column_one_based = EXTCAT_RA_COLUMN_ONE_BASED;
    std::size_t extcat_dec_column_one_based = EXTCAT_DEC_COLUMN_ONE_BASED;
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
