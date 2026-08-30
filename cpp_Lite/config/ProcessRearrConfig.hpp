#ifndef PROCESS_REARR_CONFIG_HPP
#define PROCESS_REARR_CONFIG_HPP

#include "ProcessConfig.hpp"
#include "ExtCatConfig.hpp"
#include "pathconfig.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace ProcessRearrConfig {

// ==========================================
// Configuration: Derived _all.cat column layout
// Method: Treat ichi2 as the 25 fields appended by process_main after the
//         fixed EXPO_NUM and CCD_NUM prefix columns.
// ==========================================
inline constexpr std::size_t ichi2 =
    static_cast<std::size_t>(LensingConfig::ichi2) + 1;  // Appended process_main field count.
inline constexpr std::size_t EXPO_COLUMN_COUNT = 1;  // EXPO_NUM field count.
inline constexpr std::size_t CCD_COLUMN_COUNT = 1;  // CCD_NUM field count.
inline constexpr std::size_t ALL_CAT_TOTAL_COLUMNS =
    ExtCatConfig::EXTCAT_TOTAL_COLUMNS + EXPO_COLUMN_COUNT
    + CCD_COLUMN_COUNT + ichi2;  // Default complete row width.

// ==========================================
// Configuration: Spatial partitioning and output defaults
// Method: Preserve the F77 0.1-degree full-sky grid and approximately
//         500,000 rows per deterministic weighted k-d partition.
// ==========================================
inline constexpr double SKY_GRID_DEGREES = 0.1;  // Full-sky tile width in degrees.
inline constexpr int RA_BIN_COUNT = 3600;  // Number of right-ascension bins.
inline constexpr int DEC_BIN_COUNT = 1800;  // Number of declination bins.
inline constexpr std::size_t SKY_TILE_COUNT =
    static_cast<std::size_t>(RA_BIN_COUNT) * DEC_BIN_COUNT;  // Total full-sky tile count.
inline constexpr std::uint64_t TARGET_SUBCAT_ROWS = 500000;  // Target rows per partition.
inline constexpr int SUBCAT_ID_WIDTH = 6;  // Minimum zero-padded partition ID width.
inline constexpr int OUTPUT_PRECISION = 10;  // Significant digits in catalog rows.
inline constexpr int SUMMARY_PRECISION = 4;  // Decimal places in summary bounds.
inline constexpr bool SKIP_MISSING_CATALOGS = true;  // Continue past absent input catalogs.
inline constexpr bool SKIP_MALFORMED_ROWS = true;  // Continue past malformed catalog rows.

// ==========================================
// Function: Determine the external-field width emitted by process_extcat
// Method: Use the configured pass-through width, or the explicit projection
//         length because projected output contains exactly those fields.
// ==========================================
inline std::size_t externalCatalogColumns(
    const ProcessConfig::RuntimeOptions& options) {
    return options.catalog.use_explicit_columns
               ? options.catalog.input_columns_one_based.size()
               : ExtCatConfig::EXTCAT_TOTAL_COLUMNS;
}

// ==========================================
// Function: Compute the complete _all.cat row width
// Method: Add the fixed exposure and CCD columns after the runtime-effective
//         external catalog prefix, followed by the process_main payload.
// ==========================================
inline std::size_t allCatalogColumns(
    const ProcessConfig::RuntimeOptions& options) {
    return externalCatalogColumns(options) + EXPO_COLUMN_COUNT
           + CCD_COLUMN_COUNT + ichi2;
}

static_assert(ichi2 == 25, "process_main must append 25 fields through exposure Chi2");
static_assert(ExtCatConfig::EXTCAT_TOTAL_COLUMNS > 0,
              "the external catalog column count must be positive");
static_assert(SKY_TILE_COUNT
                  <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "the full-sky MPI reduction count must fit int");

}  // namespace ProcessRearrConfig

#endif  // PROCESS_REARR_CONFIG_HPP
