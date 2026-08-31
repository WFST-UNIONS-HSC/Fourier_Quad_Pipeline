#ifndef PROCESS_CONFIG_HPP
#define PROCESS_CONFIG_HPP

#include "AstroCatConfig.hpp"
#include "ExtCatConfig.hpp"
#include "InitConfig.hpp"
#include "pathconfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ProcessConfig {

inline constexpr bool RUN_PROCESS_ASTROCAT = false;  // Run Gaia-catalog tiling by default.
inline constexpr bool RUN_PROCESS_EXTCAT = false;  // Run external-catalog tiling by default.
inline constexpr bool RUN_PROCESS_INIT = true;     // Run archive initialization by default.
inline constexpr bool RUN_PROCESS_MAIN = true;     // Run the nine-stage pipeline by default.
inline constexpr bool RUN_PROCESS_REARR = true;    // Run catalog rearrangement by default.
inline constexpr bool RUN_PROCESS_FD = false;       // Run the field-distortion test by default.

// ==========================================
// Configuration: Runtime workflow options
// Method: Seed every optional command-line value from the per-process config defaults, then
//         let the unified parser override only explicitly supplied options.
// ==========================================
struct WorkflowOptions {
    bool run_astrocat = RUN_PROCESS_ASTROCAT;  // Enable process_astrocat for this run.
    bool run_extcat = RUN_PROCESS_EXTCAT;  // Enable process_extcat for this run.
    bool run_init = RUN_PROCESS_INIT;      // Enable process_init for this run.
    bool run_main = RUN_PROCESS_MAIN;      // Enable process_main for this run.
    bool run_rearr = RUN_PROCESS_REARR;    // Enable process_rearr for this run.
    bool run_fd = RUN_PROCESS_FD;          // Enable process_fd for this run.
    bool help_requested = false;           // Print usage without running phases.
};

// ==========================================
// Configuration: Runtime process_astrocat options
// Method: Seed the independent raw-input and tile-output contract from the
//         dedicated config without coupling it to process_main at runtime.
// ==========================================
struct AstroCatOptions {
    std::string input_directory =
        AstroCatConfig::ASTROCAT_INPUT_DIRECTORY;  // Raw two-column Gaia files.
    std::string output_directory =
        AstroCatConfig::ASTROCAT_OUTPUT_DIRECTORY;  // Generated one-degree tiles.
    bool add_header = AstroCatConfig::ASTROCAT_ADD_HEADER;  // Raw input starts with data.
    std::string existing_policy =
        AstroCatConfig::ASTROCAT_EXISTING_POLICY;  // fail or overwrite.
};

struct PipelineOptions {
    std::string output_root = InitConfig::OUTPUT_ROOT;  // Shared dataset output root.
    std::vector<InitConfig::DatasetSpec> datasets = InitConfig::DATASETS;  // Sequential datasets.
    std::string exposure_list = EXPO_LIST;  // Downstream top-level exposure list.
    bool external_exposure_list_supplied = false;  // Track an explicit list override.
};

struct CatalogOptions {
    std::string directory = ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY;  // Tile output and main input.
    bool use_explicit_columns = ExtCatConfig::EXTCAT_USE_EXPLICIT_COLUMNS;  // Enable projection.
    std::vector<std::size_t> input_columns_one_based =
        ExtCatConfig::EXTCAT_INPUT_COLUMNS_ONE_BASED;  // Ordered raw-column projection.
    bool use_explicit_coordinate_columns =
        ExtCatConfig::EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS;  // Override RA/Dec discovery.
    std::size_t ra_column_one_based = ExtCatConfig::EXTCAT_RA_COLUMN_ONE_BASED;  // Raw RA field.
    std::size_t dec_column_one_based = ExtCatConfig::EXTCAT_DEC_COLUMN_ONE_BASED;  // Raw Dec field.
    std::size_t zp_column_one_based = ExtCatConfig::EXTCAT_ZP_COLUMN_ONE_BASED;  // Raw photo-z field.
};

struct ExtCatOptions {
    std::string input_directory = ExtCatConfig::EXTCAT_INPUT_DIRECTORY;  // Raw catalog root.
    std::vector<std::string> filename_tokens = ExtCatConfig::EXTCAT_FILENAME_TOKENS;  // OR filters.
    bool recursive = ExtCatConfig::EXTCAT_RECURSIVE;  // Recurse below the raw catalog root.
    std::string delimiter = ExtCatConfig::EXTCAT_DELIMITER;  // Input delimiter mode.
    std::string header_mode = ExtCatConfig::EXTCAT_HEADER_MODE;  // Header handling mode.
    std::string malformed_policy = ExtCatConfig::EXTCAT_MALFORMED_POLICY;  // Bad-row policy.
    std::string existing_policy = ExtCatConfig::EXTCAT_EXISTING_POLICY;  // Existing-tile policy.
    std::uint64_t chunk_mib = ExtCatConfig::EXTCAT_CHUNK_MIB;  // MPI task size in MiB.
};

struct InitOptions {
    std::string science_root = InitConfig::SCIENCE_ROOT;  // Science archive root.
    std::string dq_root = InitConfig::DQ_ROOT;  // DQ-mask archive root.
    std::vector<std::string> contains = InitConfig::CONTAINS;  // Archive basename filters.
    std::string existing = InitConfig::EXISTING;  // Existing-output policy.
    int f77_max_path = InitConfig::F77_MAX_PATH;  // Optional legacy path-length guard.
};

struct RearrOptions {
    std::string output_directory = REARR_OUTPUT_DIRECTORY;  // Rearranged catalog directory.
    std::string output_base_directory = REARR_OUTPUT_BASE_DIRECTORY;  // Optional output base.
    std::string exposure_list_filename = REARRANGED_EXPO_LIST_FILENAME;  // Published list name.
    std::string exposure_list_directory = REARRANGED_EXPO_LIST_DIRECTORY;  // Published list directory.
};

struct FDOptions {
    std::string exposure_list = FD_EXPO_LIST;  // Optional FD-specific exposure list.
    std::string output_directory = FD_OUTPUT_DIRECTORY;  // FD result directory.
    std::string output_base_directory = FD_OUTPUT_BASE_DIRECTORY;  // Optional FD output base.
};

struct RuntimeOptions {
    WorkflowOptions workflow;  // Phase switches and help state.
    AstroCatOptions astrocat;  // Independent Gaia-catalog tiling options.
    PipelineOptions pipeline;  // Dataset-level shared inputs and outputs.
    CatalogOptions catalog;    // Shared external-catalog contract.
    ExtCatOptions extcat;      // process_extcat parsing and discovery options.
    InitOptions init;          // process_init archive options.
    RearrOptions rearr;        // process_rearr path options.
    FDOptions fd;              // process_fd path options.
};

}  // namespace ProcessConfig

#endif  // PROCESS_CONFIG_HPP
