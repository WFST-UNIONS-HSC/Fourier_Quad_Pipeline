#ifndef EXTCAT_CONFIG_HPP
#define EXTCAT_CONFIG_HPP

// ==========================================
// ExtCatConfig - External source-catalog repartitioning defaults
// Method: Derive the tile output from the main pipeline catalog default while keeping raw-input,
//         parsing, and optional projection controls together for the first workflow phase.
//
// Note: When EXTCAT_USE_EXPLICIT_COLUMNS is true, process_main (via ExternalCatalogReader)
//       auto-indexes RA, Dec, and ZP from the ordered EXTCAT_INPUT_COLUMNS_ONE_BASED
//       projection instead of using EXTCAT_RA_COLUMN_ONE_BASED, EXTCAT_DEC_COLUMN_ONE_BASED,
//       and EXTCAT_ZP_COLUMN_ONE_BASED directly. The three field columns must still appear
//       in the projection list; their output positions follow the projection order.
// ==========================================

#include "LensingConfig.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ExtCatConfig {

inline constexpr const char* EXTCAT_INPUT_DIRECTORY = "";  // Root containing raw catalog files.
inline constexpr const char* EXTCAT_OUTPUT_DIRECTORY =
    LensingConfig::SOURCE_CAT_DEFAULT;  // Tile output and process_main input directory.
inline const std::vector<std::string> EXTCAT_FILENAME_TOKENS = {};  // OR-matched basename filters.
inline constexpr bool EXTCAT_RECURSIVE = true;  // Recurse below the input directory.
inline constexpr const char* EXTCAT_DELIMITER = "auto";  // Input delimiter detection mode.
inline constexpr const char* EXTCAT_HEADER_MODE = "auto";  // Input header handling mode.
inline constexpr const char* EXTCAT_MALFORMED_POLICY = "fail";  // Malformed-row handling policy.
inline constexpr const char* EXTCAT_EXISTING_POLICY = "fail";  // Existing-tile handling policy.
inline constexpr std::uint64_t EXTCAT_CHUNK_MIB = 64;  // MPI byte-range task size in MiB.
inline constexpr std::size_t EXTCAT_TOTAL_COLUMNS = 18;  // Canonical external catalog width.
inline constexpr bool EXTCAT_USE_EXPLICIT_COLUMNS = false;  // Enable ordered column projection.
inline const std::vector<std::size_t> EXTCAT_INPUT_COLUMNS_ONE_BASED = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
};  // Raw one-based columns emitted in output order.
inline constexpr bool EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS = false;  // Override RA/Dec discovery.
inline constexpr std::size_t EXTCAT_RA_COLUMN_ONE_BASED = 5;  // Raw one-based RA column.
inline constexpr std::size_t EXTCAT_DEC_COLUMN_ONE_BASED = 6;  // Raw one-based Dec column.
inline constexpr std::size_t EXTCAT_ZP_COLUMN_ONE_BASED = 17;  // Raw one-based photo-z column.

}  // namespace ExtCatConfig

#endif  // EXTCAT_CONFIG_HPP
