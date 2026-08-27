#include "general/MPIScheduler.hpp"
#include "general/NumericalRecipes.hpp"
#include "general/ExposureList.hpp"
#include "general/PathUtils.hpp"
#include "process_main/ExternalCatalogReader.hpp"
#include "LensingConfig.hpp"
#include "ProcessConfig.hpp"
#include "ExtCatConfig.hpp"
#include "InitConfig.hpp"
#include "process_extcat/process_extcat.hpp"
#include "process_init/process_init.hpp"
#include "process_main/process_main.hpp"
#include "process_rearr/CatalogRearranger.hpp"
#include "process_rearr/process_rearr.hpp"
#include "process_fd/process_fd.hpp"
#include "FDConfig.hpp"
#include "ProcessRearrConfig.hpp"

#include <mpi.h>

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParserState {
    bool extcat_contains_option_seen = false;
    bool dataset_option_seen = false;
    bool legacy_dataset_option_seen = false;
    bool contains_option_seen = false;
};

// ==========================================
// Function: Parse one command-line boolean
// Method: Accept common textual and numeric forms without locale-dependent conversion.
// ==========================================
bool parseBoolean(const std::string& value, bool& parsed) {
    if (value == "true" || value == "1" || value == "on") {
        parsed = true;
        return true;
    }
    if (value == "false" || value == "0" || value == "off") {
        parsed = false;
        return true;
    }
    return false;
}

// ==========================================
// Function: Parse one non-negative command-line integer
// Method: Require full decimal consumption and guard the C++ int range.
// ==========================================
bool parseNonNegativeInteger(const std::string& value, int& parsed) {
    std::size_t consumed = 0;
    long long number = 0;
    try {
        number = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number < 0
        || number > std::numeric_limits<int>::max()) {
        return false;
    }
    parsed = static_cast<int>(number);
    return true;
}

// ==========================================
// Function: Parse one positive unsigned command-line integer
// Method: Require full decimal consumption and guard the uint64 range.
// ==========================================
bool parsePositiveUnsigned(const std::string& value, std::uint64_t& parsed) {
    if (value.empty() || value.front() == '-') {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number == 0
        || number > std::numeric_limits<std::uint64_t>::max()) {
        return false;
    }
    parsed = static_cast<std::uint64_t>(number);
    return true;
}

// ==========================================
// Function: Parse the external-catalog input column projection
// Method: Read one or more comma-separated positive one-based indices in output order.
// ==========================================
bool parseExtcatColumns(const std::string& value,
                        std::vector<std::size_t>& columns) {
    std::stringstream stream(value);
    std::string token;
    std::vector<std::size_t> parsed;
    while (std::getline(stream, token, ',')) {
        if (token.empty() || token.front() == '-') {
            return false;
        }
        std::size_t consumed = 0;
        unsigned long long number = 0;
        try {
            number = std::stoull(token, &consumed);
        } catch (const std::exception&) {
            return false;
        }
        if (consumed != token.size() || number == 0
            || number > std::numeric_limits<std::size_t>::max()) {
            return false;
        }
        parsed.push_back(static_cast<std::size_t>(number));
    }
    if (parsed.empty()) {
        return false;
    }
    columns = parsed;
    return true;
}

// ==========================================
// Function: Parse one external-catalog field column
// Method: Require one positive one-based RA, Dec, or ZP index within the platform size_t range.
// ==========================================
bool parseExtcatFieldColumn(const std::string& value, std::size_t& column) {
    if (value.empty() || value.front() == '-') {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number == 0
        || number > std::numeric_limits<std::size_t>::max()) {
        return false;
    }
    column = static_cast<std::size_t>(number);
    return true;
}

// ==========================================
// Function: Parse one paired dataset value
// Method: Split exactly one TARGET:PREFIX value and reject missing components.
// ==========================================
bool parseDataset(const std::string& value,
                  InitConfig::DatasetSpec& dataset,
                  std::string& error) {
    const std::size_t separator = value.find(':');
    if (separator == std::string::npos || separator == 0 || separator + 1 == value.size()
        || value.find(':', separator + 1) != std::string::npos) {
        error = "--dataset must use TARGET:PREFIX with both components non-empty";
        return false;
    }
    dataset.target = value.substr(0, separator);
    dataset.prefix = value.substr(separator + 1);
    return true;
}

// ==========================================
// Function: Enter legacy single-dataset override mode
// Method: Preserve the first configured pair, collapse the list to one entry,
//         and forbid mixing legacy target/prefix flags with --dataset.
// ==========================================
bool prepareLegacyDataset(ProcessConfig::RuntimeOptions& options,
                          ParserState& state,
                          std::string& error) {
    if (state.dataset_option_seen) {
        error = "--dataset cannot be combined with --target or --prefix";
        return false;
    }
    if (!state.legacy_dataset_option_seen) {
        InitConfig::DatasetSpec dataset;
        if (!options.pipeline.datasets.empty()) {
            dataset = options.pipeline.datasets.front();
        }
        options.pipeline.datasets.assign(1, dataset);
        state.legacy_dataset_option_seen = true;
    }
    return true;
}

// ==========================================
// Function: Apply one normalized named option
// Method: Override scalar defaults, accumulate batch values, and validate typed
//         values immediately while tracking incompatible dataset syntaxes.
// ==========================================
bool applyNamedOption(const std::string& name,
                      const std::string& value,
                      ProcessConfig::RuntimeOptions& options,
                      ParserState& state,
                      std::string& error) {
    if (name == "--run-extcat") {
        if (!parseBoolean(value, options.workflow.run_extcat)) {
            error = "--run-extcat must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--extcat-input") {
        options.extcat.input_directory = value;
    } else if (name == "--extcat-output") {
        options.catalog.directory = value;
    } else if (name == "--extcat-contains") {
        if (value.empty()) {
            error = "--extcat-contains must not be empty";
            return false;
        }
        if (!state.extcat_contains_option_seen) {
            options.extcat.filename_tokens.clear();
            state.extcat_contains_option_seen = true;
        }
        options.extcat.filename_tokens.push_back(value);
    } else if (name == "--extcat-recursive") {
        if (!parseBoolean(value, options.extcat.recursive)) {
            error = "--extcat-recursive must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--extcat-delimiter") {
        if (value != "auto" && value != "whitespace" && value != "comma"
            && value != "tab") {
            error = "--extcat-delimiter must be auto, whitespace, comma, or tab";
            return false;
        }
        options.extcat.delimiter = value;
    } else if (name == "--extcat-header") {
        if (value != "auto" && value != "present" && value != "absent") {
            error = "--extcat-header must be auto, present, or absent";
            return false;
        }
        options.extcat.header_mode = value;
    } else if (name == "--extcat-columns") {
        if (!parseExtcatColumns(value, options.catalog.input_columns_one_based)) {
            error = "--extcat-columns must contain one or more positive one-based indices";
            return false;
        }
        options.catalog.use_explicit_columns = true;
    } else if (name == "--extcat-ra-column") {
        if (!parseExtcatFieldColumn(value, options.catalog.ra_column_one_based)) {
            error = "--extcat-ra-column must be a positive one-based index";
            return false;
        }
        options.catalog.use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-dec-column") {
        if (!parseExtcatFieldColumn(value, options.catalog.dec_column_one_based)) {
            error = "--extcat-dec-column must be a positive one-based index";
            return false;
        }
        options.catalog.use_explicit_coordinate_columns = true;
    } else if (name == "--extcat-zp-column") {
        if (!parseExtcatFieldColumn(value, options.catalog.zp_column_one_based)) {
            error = "--extcat-zp-column must be a positive one-based index";
            return false;
        }
    } else if (name == "--extcat-chunk-mib") {
        if (!parsePositiveUnsigned(value, options.extcat.chunk_mib)) {
            error = "--extcat-chunk-mib must be a positive integer";
            return false;
        }
    } else if (name == "--extcat-malformed") {
        if (value != "fail" && value != "skip") {
            error = "--extcat-malformed must be fail or skip";
            return false;
        }
        options.extcat.malformed_policy = value;
    } else if (name == "--extcat-existing") {
        if (value != "fail" && value != "overwrite") {
            error = "--extcat-existing must be fail or overwrite";
            return false;
        }
        options.extcat.existing_policy = value;
    } else if (name == "--run-init") {
        if (!parseBoolean(value, options.workflow.run_init)) {
            error = "--run-init must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-main") {
        if (!parseBoolean(value, options.workflow.run_main)) {
            error = "--run-main must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-rearr") {
        if (!parseBoolean(value, options.workflow.run_rearr)) {
            error = "--run-rearr must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-fd") {
        if (!parseBoolean(value, options.workflow.run_fd)) {
            error = "--run-fd must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--science-root") {
        options.init.science_root = value;
    } else if (name == "--dq-root") {
        options.init.dq_root = value;
    } else if (name == "--output-root") {
        options.pipeline.output_root = value;
    } else if (name == "--dataset") {
        if (state.legacy_dataset_option_seen) {
            error = "--dataset cannot be combined with --target or --prefix";
            return false;
        }
        InitConfig::DatasetSpec dataset;
        if (!parseDataset(value, dataset, error)) {
            return false;
        }
        if (!state.dataset_option_seen) {
            options.pipeline.datasets.clear();
            state.dataset_option_seen = true;
        }
        options.pipeline.datasets.push_back(dataset);
    } else if (name == "--target") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.pipeline.datasets.front().target = value;
    } else if (name == "--prefix") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.pipeline.datasets.front().prefix = value;
    } else if (name == "--contains") {
        if (value.empty()) {
            error = "--contains must not be empty";
            return false;
        }
        if (!state.contains_option_seen) {
            options.init.contains.clear();
            state.contains_option_seen = true;
        }
        options.init.contains.push_back(value);
    } else if (name == "--existing") {
        if (value != "fail" && value != "resume" && value != "overwrite") {
            error = "--existing must be fail, resume, or overwrite";
            return false;
        }
        options.init.existing = value;
    } else if (name == "--f77-max-path") {
        if (!parseNonNegativeInteger(value, options.init.f77_max_path)) {
            error = "--f77-max-path must be a non-negative integer";
            return false;
        }
    } else if (name == "--expo-list") {
        options.pipeline.exposure_list = value;
        options.pipeline.external_exposure_list_supplied = true;
    } else if (name == "--rearr-output-dir") {
        options.rearr.output_directory = value;
    } else if (name == "--rearr-output-base") {
        options.rearr.output_base_directory = value;
    } else if (name == "--rearr-list-name") {
        options.rearr.exposure_list_filename = value;
    } else if (name == "--rearr-list-dir") {
        options.rearr.exposure_list_directory = value;
    } else if (name == "--fd-expo-list") {
        options.fd.exposure_list = value;
    } else if (name == "--fd-output-dir") {
        options.fd.output_directory = value;
    } else if (name == "--fd-output-base") {
        options.fd.output_base_directory = value;
    } else {
        error = "unknown option: " + name;
        return false;
    }
    return true;
}

// ==========================================
// Function: Validate the effective workflow configuration
// Method: Enforce executable modes, unique safe targets, paired prefixes, valid
//         token lists, and unambiguous exposure-list input for downstream-only runs.
// ==========================================
bool validateOptions(const ProcessConfig::RuntimeOptions& options, std::string& error) {
    if (options.workflow.help_requested) {
        return true;
    }
    if (!options.workflow.run_extcat && !options.workflow.run_init
        && !options.workflow.run_fd
        && !options.workflow.run_main && !options.workflow.run_rearr) {
        error = "--run-extcat, --run-init, --run-main, --run-rearr, and --run-fd cannot all be false";
        return false;
    }
    if ((options.workflow.run_extcat || options.workflow.run_main)
        && options.catalog.directory.empty()) {
        error = "external source-catalog output directory must not be empty";
        return false;
    }
    if (options.workflow.run_extcat && options.extcat.input_directory.empty()) {
        error = "external source-catalog input directory must not be empty";
        return false;
    }
    if ((options.workflow.run_init || options.workflow.run_main || options.workflow.run_fd
         || options.workflow.run_rearr)
        && options.pipeline.datasets.empty()) {
        error = "at least one dataset must be configured or supplied with --dataset";
        return false;
    }

    std::set<std::string> targets;
    if (options.workflow.run_init || options.workflow.run_main || options.workflow.run_fd
        || options.workflow.run_rearr) {
        for (const InitConfig::DatasetSpec& dataset : options.pipeline.datasets) {
            if (dataset.target.empty() || dataset.target == "." || dataset.target == ".."
                || dataset.target.find('/') != std::string::npos
                || dataset.target.find('\\') != std::string::npos) {
                error = "each dataset target must be one non-empty directory name";
                return false;
            }
            if (dataset.prefix.empty()) {
                error = "each dataset prefix must be non-empty";
                return false;
            }
            if (!targets.insert(dataset.target).second) {
                error = "dataset target is duplicated: " + dataset.target;
                return false;
            }
        }
    }

    if (options.workflow.run_init) {
        for (const std::string& token : options.init.contains) {
            if (token.empty()) {
                error = "contains tokens must be non-empty";
                return false;
            }
        }
    }
    if (options.workflow.run_extcat) {
        for (const std::string& token : options.extcat.filename_tokens) {
            if (token.empty()) {
                error = "external-catalog contains tokens must be non-empty";
                return false;
            }
        }
        if (options.catalog.use_explicit_columns
            && options.catalog.input_columns_one_based.empty()) {
            error = "external-catalog explicit column list must not be empty";
            return false;
        }
        if (options.catalog.use_explicit_coordinate_columns
            && (options.catalog.ra_column_one_based == 0
                || options.catalog.dec_column_one_based == 0
                || options.catalog.ra_column_one_based
                       == options.catalog.dec_column_one_based)) {
            error = "external-catalog RA and Dec columns must be distinct positive indices";
            return false;
        }
    }
    if (options.workflow.run_main) {
        ExternalCatalogReader::ColumnSelection selection;
        if (!ExternalCatalogReader::resolveColumnSelection(options, selection, error)) {
            return false;
        }
    }
    if (options.workflow.run_rearr) {
        ProcessRearr::CatalogLayout layout;
        if (!ProcessRearr::resolveCatalogLayout(options, layout, error)) {
            return false;
        }
    }

    if ((options.workflow.run_main || options.workflow.run_rearr)
        && !options.workflow.run_init
        && options.pipeline.datasets.size() > 1 && !options.pipeline.exposure_list.empty()) {
        error = "one --expo-list cannot serve multiple datasets in downstream-only mode; "
                "omit it to derive expo_<target>.list for each dataset";
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the unified variable-length workflow command line
// Method: Accept --name value and --name=value in any order, accumulate dataset
//         and token options, and retain one positional exposure-list alias.
// ==========================================
bool parseCommandLine(int argc,
                      char** argv,
                      ProcessConfig::RuntimeOptions& options,
                      std::string& error) {
    std::string legacy_exposure_list;
    bool legacy_exposure_list_supplied = false;
    ParserState state;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            options.workflow.help_requested = true;
            continue;
        }

        if (argument.rfind("--", 0) != 0) {
            if (legacy_exposure_list_supplied) {
                error = "only one positional exposure-list compatibility argument is allowed";
                return false;
            }
            legacy_exposure_list = argument;
            legacy_exposure_list_supplied = true;
            continue;
        }

        const std::size_t equals = argument.find('=');
        const std::string name = argument.substr(0, equals);
        std::string value;
        if (equals != std::string::npos) {
            value = argument.substr(equals + 1);
        } else {
            if (index + 1 >= argc || std::string(argv[index + 1]).rfind("--", 0) == 0) {
                error = "missing value after " + name;
                return false;
            }
            value = argv[++index];
        }

        if (!applyNamedOption(name, value, options, state, error)) {
            return false;
        }
    }

    if (!options.pipeline.external_exposure_list_supplied && legacy_exposure_list_supplied) {
        options.pipeline.exposure_list = legacy_exposure_list;
        options.pipeline.external_exposure_list_supplied = true;
    }
    return validateOptions(options, error);
}

// ==========================================
// Function: Resolve one dataset exposure list used by downstream-only mode
// Method: Prefer the single external/default list, otherwise derive the current
//         dataset's expo_<target>.list and normalize it to an absolute path.
// ==========================================
std::string resolveExposureList(const ProcessConfig::RuntimeOptions& options,
                                const InitConfig::DatasetSpec& dataset) {
    std::filesystem::path path;
    if (!options.pipeline.exposure_list.empty()) {
        path = options.pipeline.exposure_list;
    } else {
        if (options.pipeline.output_root.empty()) {
            throw std::runtime_error(
                "--expo-list is absent and output-root cannot derive its default");
        }
        path = std::filesystem::path(options.pipeline.output_root)
               / ("expo_" + dataset.target + ".list");
    }
    return PathUtils::normalizedAbsolute(path).string();
}

// ==========================================
// Function: Format configured dataset defaults for help output
// Method: Join every target/prefix pair without mutating the configured list.
// ==========================================
std::string configuredDatasetsText() {
    if (InitConfig::DATASETS.empty()) {
        return "none";
    }
    std::string text;
    for (const InitConfig::DatasetSpec& dataset : InitConfig::DATASETS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += dataset.target + ":" + dataset.prefix;
    }
    return text;
}

// ==========================================
// Function: Format configured token defaults for help output
// Method: Join every OR-matched basename token, or state that filtering is disabled.
// ==========================================
std::string configuredContainsText() {
    if (InitConfig::CONTAINS.empty()) {
        return "none (no token filter)";
    }
    std::string text;
    for (const std::string& token : InitConfig::CONTAINS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += token;
    }
    return text;
}

// ==========================================
// Function: Format external-catalog filename-token defaults
// Method: Join OR-matched raw-catalog basename tokens or state that all files match.
// ==========================================
std::string configuredExtcatContainsText() {
    if (ExtCatConfig::EXTCAT_FILENAME_TOKENS.empty()) {
        return "none (all files)";
    }
    std::string text;
    for (const std::string& token : ExtCatConfig::EXTCAT_FILENAME_TOKENS) {
        if (!text.empty()) {
            text += ", ";
        }
        text += token;
    }
    return text;
}

// ==========================================
// Function: Print the unified workflow command-line contract
// Method: Describe batch/list accumulation, initializer values, list precedence,
//         compatibility input, and the configured defaults.
// ==========================================
void printUsage(const char* program_name) {
    std::cout
        << "Usage: " << program_name << " [options] [LEGACY_EXPO_LIST]\n"
        << "  --run-extcat BOOL     Repartition raw external catalogs first (default: "
        << (ProcessConfig::RUN_PROCESS_EXTCAT ? "true" : "false") << ")\n"
        << "  --run-init BOOL       Run initializer (default: "
        << (ProcessConfig::RUN_PROCESS_INIT ? "true" : "false") << ")\n"
        << "  --run-main BOOL       Run numerical pipeline (default: "
        << (ProcessConfig::RUN_PROCESS_MAIN ? "true" : "false") << ")\n"
        << "  --run-rearr BOOL      Rearrange _all.cat; follows process_main if both run (default: "
        << (ProcessConfig::RUN_PROCESS_REARR ? "true" : "false") << ")\n"
        << "  --run-fd BOOL         Run FD (field-distortion) shear test; follows process_main (default: "
        << (ProcessConfig::RUN_PROCESS_FD ? "true" : "false") << ")\n"
        << "  --extcat-input PATH   Directory containing raw external catalogs\n"
        << "  --extcat-output PATH  External-catalog tile directory used by process_main\n"
        << "  --extcat-contains T   Repeatable raw basename token, matched with OR (default: "
        << configuredExtcatContainsText() << ")\n"
        << "  --extcat-recursive B  Recurse below extcat input (default: "
        << (ExtCatConfig::EXTCAT_RECURSIVE ? "true" : "false") << ")\n"
        << "  --extcat-delimiter M  auto, whitespace, comma, or tab\n"
        << "  --extcat-header M     auto, present, or absent\n"
        << "  --extcat-columns LIST Ordered one-based input indices; output width follows LIST\n"
        << "  --extcat-ra-column N  Raw one-based RA column; overrides header discovery\n"
        << "  --extcat-dec-column N Raw one-based Dec column; overrides header discovery\n"
        << "  --extcat-zp-column N  Raw one-based ZP column consumed by process_main\n"
        << "  --extcat-chunk-mib N  MPI byte-range task size in MiB (default: "
        << ExtCatConfig::EXTCAT_CHUNK_MIB << ")\n"
        << "  --extcat-malformed P  fail or skip malformed rows\n"
        << "  --extcat-existing P   fail or overwrite generated tiles\n"
        << "  --science-root PATH   Original Science FITS/FZ repository\n"
        << "  --dq-root PATH        Original DQ FITS/FZ repository\n"
        << "  --output-root PATH    Parent of targets and generated exposure lists\n"
        << "  --dataset T:P         Repeatable TARGET:PREFIX dataset pair (default: "
        << configuredDatasetsText() << ")\n"
        << "  --target NAME         Legacy single-dataset target; do not mix with --dataset\n"
        << "  --prefix TEXT         Legacy single-dataset prefix; do not mix with --dataset\n"
        << "  --contains TEXT       Repeatable basename token, matched with OR (default: "
        << configuredContainsText() << ")\n"
        << "  --existing MODE       fail, resume, or overwrite (default: "
        << InitConfig::EXISTING << ")\n"
        << "  --f77-max-path N      Generated path limit; zero disables (default: "
        << InitConfig::F77_MAX_PATH << ")\n"
        << "  --expo-list PATH      Single exposure list for main/rearr-only mode\n"
        << "  --rearr-output-dir D  process_rearr output sub-directory (default: "
        << ProcessConfig::REARR_OUTPUT_DIRECTORY << ")\n"
        << "  --rearr-output-base P process_rearr output base path (default: dataset root)\n"
        << "  --rearr-list-name F   Rearranged expo-list filename (default: "
        << ProcessConfig::REARRANGED_EXPO_LIST_FILENAME << ")\n"
        << "  --rearr-list-dir P    Rearranged expo-list directory (default: expo-list parent)\n"
        << "  --fd-expo-list PATH   process_fd exposure list (default: rearranged list)\n"
        << "  --fd-output-dir D    process_fd output sub-directory (default: "
        << ProcessConfig::FD_OUTPUT_DIRECTORY << ")\n"
        << "  --fd-output-base P    process_fd output base path (default: dataset root)\n"
        << "  --help                Show this help\n"
        << "Options accept both --name value and --name=value. The first explicit "
           "--dataset, --contains, or --extcat-contains replaces its corresponding "
           "configured list; repeats append.\n"
        << "Other duplicate scalar options use the last value. Downstream-only batches derive "
           "expo_<target>.list per dataset when --expo-list is omitted.\n"
        << "Without --extcat-columns, all raw catalog fields keep their original order; "
           "otherwise output width and order follow LIST exactly.\n"
        << "When process_main runs, LIST must contain the configured RA, Dec, and ZP "
           "raw columns. A rearr-only run requires RA and Dec but not ZP; output "
           "positions follow LIST order.\n"
        << "process_extcat runs once before the dataset loop. When later phases run, each "
           "process_init-generated absolute exposure-list "
           "path overrides external input for its dataset.\n";
}


// ==========================================
// Function: Derive dataset root from the first image path in an expo list
// Method: Read the first per-exposure list entry, open it, read the first
//         non-empty image line, and compute the great-grandparent directory
//         (getDir level 3), matching process_main's dir_output derivation.
// ==========================================
std::string deriveDatasetRootFromExpoList(const std::string& exposure_list) {
    std::vector<ExposureList::Entry> entries;
    std::string exposure_error;
    if (!ExposureList::loadPipelineList(
            exposure_list, entries,
            static_cast<std::size_t>(LensingConfig::NMAX_EXPO), exposure_error)) {
        return "";
    }
    for (const ExposureList::Entry& entry : entries) {
        std::ifstream list_input(entry.path);
        if (!list_input.is_open()) {
            continue;
        }
        std::string line;
        while (std::getline(list_input, line)) {
            std::size_t first = 0;
            while (first < line.size()
                   && (line[first] == ' ' || line[first] == '\t'
                       || line[first] == '\r' || line[first] == '\n')) {
                ++first;
            }
            std::size_t last = line.size();
            while (last > first
                   && (line[last - 1] == ' ' || line[last - 1] == '\t'
                       || line[last - 1] == '\r' || line[last - 1] == '\n')) {
                --last;
            }
            if (first == last) {
                continue;
            }
            line = line.substr(first, last - first);
            if (line.size() >= 2 && line.front() == '"' && line.back() == '"') {
                line = line.substr(1, line.size() - 2);
            }
            if (line.empty()) {
                continue;
            }
            const std::filesystem::path image_path(line);
            std::filesystem::path great_grandparent;
            std::string parent_error;
            if (!PathUtils::parentAtLevel(
                    image_path, 3, great_grandparent, parent_error)) {
                return "";
            }
            return PathUtils::normalizedAbsolute(great_grandparent).string();
        }
    }
    return "";
}

// ==========================================
// Function: Parse and collectively validate runtime options
// Method: Run the existing parser on every rank, combine success with MPI_MIN,
//         and report one rank-zero usage error without exiting locally.
// ==========================================
int prepareRuntimeOptions(int argc,
                          char** argv,
                          ProcessConfig::RuntimeOptions& options) {
    std::string parse_error;
    const int local_parse_ok =
        parseCommandLine(argc, argv, options, parse_error) ? 1 : 0;
    int global_parse_ok = 0;
    MPI_Allreduce(&local_parse_ok, &global_parse_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (global_parse_ok != 0) {
        return 0;
    }
    if (MPIScheduler::state.rank == 0) {
        std::cerr << "Argument error: "
                  << (parse_error.empty()
                          ? "parsing failed on another rank"
                          : parse_error)
                  << std::endl;
        printUsage(argv[0]);
    }
    return 2;
}

// ==========================================
// Function: Run the optional external-catalog phase
// Method: Preserve its once-before-datasets position and the following
//         synchronization barrier only after success.
// ==========================================
int runExtcatPhase(const ProcessConfig::RuntimeOptions& options) {
    if (!options.workflow.run_extcat) {
        return 0;
    }
    if (MPIScheduler::state.rank == 0) {
        std::cout << "Running process_extcat before all dataset phases"
                  << std::endl;
    }
    const int result = process_extcat(options, MPI_COMM_WORLD);
    if (result == 0) {
        MPIScheduler::barrier();
    }
    return result;
}

// ==========================================
// Function: Resolve downstream exposure-list input for one dataset
// Method: Prefer successful initializer output; otherwise resolve the
//         configured or derived list collectively and preserve diagnostics.
// ==========================================
int resolveDatasetExposureList(const ProcessConfig::RuntimeOptions& options,
                               const InitConfig::DatasetSpec& dataset,
                               const std::string& generated_exposure_list,
                               std::string& selected_exposure_list) {
    if (options.workflow.run_init) {
        selected_exposure_list = generated_exposure_list;
        if (MPIScheduler::state.rank == 0) {
            if (!options.pipeline.exposure_list.empty()) {
                std::cout << "External exposure list overridden by process_init output: ";
            } else {
                std::cout << "Downstream phases will use process_init output: ";
            }
            std::cout << selected_exposure_list << std::endl;
        }
        MPIScheduler::barrier();
        return 0;
    }

    int local_path_ok = 1;
    std::string path_error;
    try {
        selected_exposure_list = resolveExposureList(options, dataset);
    } catch (const std::exception& exception) {
        local_path_ok = 0;
        path_error = exception.what();
    }
    int global_path_ok = 0;
    MPI_Allreduce(&local_path_ok, &global_path_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (global_path_ok != 0) {
        return 0;
    }
    if (MPIScheduler::state.rank == 0) {
        std::cerr << "Exposure-list argument error: "
                  << (path_error.empty()
                          ? "path resolution failed on another rank"
                          : path_error)
                  << std::endl;
    }
    return 2;
}

// ==========================================
// Function: Run process_main for one dataset
// Method: Initialize the per-rank RNG once for the workflow, preserve the
//         historical seed report/barrier, then invoke the numerical pipeline.
// ==========================================
int runMainPhase(const std::string& exposure_list,
                 const ProcessConfig::RuntimeOptions& options,
                 bool& rng_initialized) {
    if (!options.workflow.run_main) {
        return 0;
    }
    if (!rng_initialized) {
        const unsigned int rng_seed = NumericalRecipes::initializeRan1Seed(
            MPIScheduler::state.rank, MPIScheduler::state.size);
        std::cout << "RNG_SEED rank seed: " << MPIScheduler::state.rank
                  << " " << rng_seed << std::endl;
        MPIScheduler::barrier();
        rng_initialized = true;
    }
    return process_main(exposure_list, options);
}

// ==========================================
// Function: Run process_rearr for one dataset
// Method: Preserve the pre-phase barrier, rank-zero status line, and existing
//         runtime options contract.
// ==========================================
int runRearrPhase(const std::string& exposure_list,
                  const ProcessConfig::RuntimeOptions& options) {
    if (!options.workflow.run_rearr) {
        return 0;
    }
    MPIScheduler::barrier();
    if (MPIScheduler::state.rank == 0) {
        std::cout << "Running process_rearr"
                  << (options.workflow.run_main ? " after process_main" : "")
                  << std::endl;
    }
    return process_rearr(exposure_list, options, MPI_COMM_WORLD);
}

// ==========================================
// Function: Run process_fd for one dataset
// Method: Preserve list derivation, dataset-root discovery, RNG setup, and
//         phase ordering before invoking the FD measurement.
// ==========================================
int runFdPhase(const std::string& selected_exposure_list,
               const ProcessConfig::RuntimeOptions& options,
               bool& rng_initialized) {
    if (!options.workflow.run_fd) {
        return 0;
    }
    MPIScheduler::barrier();
    if (MPIScheduler::state.rank == 0) {
        std::cout << "Running process_fd" << std::endl;
    }
    if (!rng_initialized) {
        NumericalRecipes::initializeRan1Seed(
            MPIScheduler::state.rank, MPIScheduler::state.size);
        MPIScheduler::barrier();
        rng_initialized = true;
    }

    std::string fd_exposure_list;
    if (options.fd.exposure_list.empty()) {
        const std::filesystem::path list_directory =
            options.rearr.exposure_list_directory.empty()
                ? std::filesystem::path(selected_exposure_list).parent_path()
                : std::filesystem::path(options.rearr.exposure_list_directory);
        fd_exposure_list =
            std::filesystem::absolute(
                list_directory / options.rearr.exposure_list_filename)
                .lexically_normal()
                .string();
    } else {
        fd_exposure_list = options.fd.exposure_list;
    }
    const std::string dataset_root =
        deriveDatasetRootFromExpoList(selected_exposure_list);
    if (MPIScheduler::state.rank == 0) {
        std::cout << "FD expo list: " << fd_exposure_list
                  << "  dataset_root: " << dataset_root << std::endl;
    }
    return process_fd(fd_exposure_list, options, dataset_root, MPI_COMM_WORLD);
}

// ==========================================
// Function: Run one dataset's enabled phase chain
// Method: Execute init, resolve the shared exposure list, then run main,
//         rearrangement, and FD in their original fail-fast order.
// ==========================================
int runDataset(const InitConfig::DatasetSpec& dataset,
               std::size_t index,
               std::size_t dataset_count,
               const ProcessConfig::RuntimeOptions& options,
               bool& rng_initialized) {
    if (MPIScheduler::state.rank == 0) {
        std::cout << "Dataset " << (index + 1) << "/" << dataset_count
                  << ": target=" << dataset.target
                  << " prefix=" << dataset.prefix << std::endl;
    }

    std::string generated_exposure_list;
    int result = 0;
    if (options.workflow.run_init) {
        result = process_init(options, dataset, generated_exposure_list);
    }

    std::string selected_exposure_list;
    const bool downstream_enabled = options.workflow.run_main
                                    || options.workflow.run_rearr
                                    || options.workflow.run_fd;
    if (result == 0 && downstream_enabled) {
        result = resolveDatasetExposureList(
            options, dataset, generated_exposure_list, selected_exposure_list);
    }
    if (result == 0) {
        result = runMainPhase(selected_exposure_list, options, rng_initialized);
    }
    if (result == 0) {
        result = runRearrPhase(selected_exposure_list, options);
    }
    if (result == 0) {
        result = runFdPhase(selected_exposure_list, options, rng_initialized);
    }
    return result;
}

// ==========================================
// Function: Run the complete configured workflow
// Method: Execute process_extcat once, then process datasets sequentially with
//         fail-fast return codes and the existing inter-dataset barrier.
// ==========================================
int runWorkflow(const ProcessConfig::RuntimeOptions& options) {
    int result = runExtcatPhase(options);
    bool rng_initialized = false;
    const bool dataset_phase_enabled = options.workflow.run_init
                                       || options.workflow.run_main
                                       || options.workflow.run_rearr
                                       || options.workflow.run_fd;
    for (std::size_t index = 0;
         index < options.pipeline.datasets.size() && result == 0
             && dataset_phase_enabled;
         ++index) {
        result = runDataset(options.pipeline.datasets[index], index,
                            options.pipeline.datasets.size(), options,
                            rng_initialized);
        if (result == 0 && index + 1 < options.pipeline.datasets.size()) {
            MPIScheduler::barrier();
        }
    }
    return result;
}

}  // namespace

// ==========================================
// Function: Initialize, run, and finalize the unified pipeline workflow
// Method: Keep MPI ownership in one entry point and delegate parsing and phase
//         orchestration to the local helpers above.
// ==========================================
int main(int argc, char* argv[]) {
    MPIScheduler::init(argc, argv);
    if (MPIScheduler::state.rank == 0) {
        std::cout << "MPI Init Done..." << std::endl;
    }

    ProcessConfig::RuntimeOptions options;
    int result = prepareRuntimeOptions(argc, argv, options);
    if (result == 0 && options.workflow.help_requested) {
        if (MPIScheduler::state.rank == 0) {
            printUsage(argv[0]);
        }
    } else if (result == 0) {
        result = runWorkflow(options);
    }

    MPIScheduler::finalize();
    return result;
}
