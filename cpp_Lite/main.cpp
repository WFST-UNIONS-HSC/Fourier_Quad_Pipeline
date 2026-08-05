#include "process_main/MPIScheduler.hpp"
#include "process_main/NumericalRecipes.hpp"
#include "ProcessConfig.hpp"
#include "process_init/process_init.hpp"
#include "process_main/process_main.hpp"

#include <mpi.h>

#include <filesystem>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParserState {
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
// Function: Parse one paired dataset value
// Method: Split exactly one TARGET:PREFIX value and reject missing components.
// ==========================================
bool parseDataset(const std::string& value,
                  ProcessConfig::DatasetSpec& dataset,
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
        ProcessConfig::DatasetSpec dataset;
        if (!options.datasets.empty()) {
            dataset = options.datasets.front();
        }
        options.datasets.assign(1, dataset);
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
    if (name == "--run-init") {
        if (!parseBoolean(value, options.run_process_init)) {
            error = "--run-init must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--run-main") {
        if (!parseBoolean(value, options.run_process_main)) {
            error = "--run-main must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--science-root") {
        options.science_root = value;
    } else if (name == "--dq-root") {
        options.dq_root = value;
    } else if (name == "--output-root") {
        options.output_root = value;
    } else if (name == "--dataset") {
        if (state.legacy_dataset_option_seen) {
            error = "--dataset cannot be combined with --target or --prefix";
            return false;
        }
        ProcessConfig::DatasetSpec dataset;
        if (!parseDataset(value, dataset, error)) {
            return false;
        }
        if (!state.dataset_option_seen) {
            options.datasets.clear();
            state.dataset_option_seen = true;
        }
        options.datasets.push_back(dataset);
    } else if (name == "--target") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.datasets.front().target = value;
    } else if (name == "--prefix") {
        if (!prepareLegacyDataset(options, state, error)) {
            return false;
        }
        options.datasets.front().prefix = value;
    } else if (name == "--contains") {
        if (value.empty()) {
            error = "--contains must not be empty";
            return false;
        }
        if (!state.contains_option_seen) {
            options.contains.clear();
            state.contains_option_seen = true;
        }
        options.contains.push_back(value);
    } else if (name == "--existing") {
        if (value != "fail" && value != "resume" && value != "overwrite") {
            error = "--existing must be fail, resume, or overwrite";
            return false;
        }
        options.existing = value;
    } else if (name == "--f77-max-path") {
        if (!parseNonNegativeInteger(value, options.f77_max_path)) {
            error = "--f77-max-path must be a non-negative integer";
            return false;
        }
    } else if (name == "--expo-list") {
        options.expo_list = value;
        options.external_expo_list_supplied = true;
    } else {
        error = "unknown option: " + name;
        return false;
    }
    return true;
}

// ==========================================
// Function: Validate the effective workflow configuration
// Method: Enforce executable modes, unique safe targets, paired prefixes, valid
//         token lists, and unambiguous exposure-list input for batch main-only runs.
// ==========================================
bool validateOptions(const ProcessConfig::RuntimeOptions& options, std::string& error) {
    if (options.help_requested) {
        return true;
    }
    if (!options.run_process_init && !options.run_process_main) {
        error = "--run-init and --run-main cannot both be false";
        return false;
    }
    if (options.datasets.empty()) {
        error = "at least one dataset must be configured or supplied with --dataset";
        return false;
    }

    std::set<std::string> targets;
    for (const ProcessConfig::DatasetSpec& dataset : options.datasets) {
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

    for (const std::string& token : options.contains) {
        if (token.empty()) {
            error = "contains tokens must be non-empty";
            return false;
        }
    }

    if (options.run_process_main && !options.run_process_init
        && options.datasets.size() > 1 && !options.expo_list.empty()) {
        error = "one --expo-list cannot serve multiple datasets in main-only mode; "
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
            options.help_requested = true;
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

    if (!options.external_expo_list_supplied && legacy_exposure_list_supplied) {
        options.expo_list = legacy_exposure_list;
        options.external_expo_list_supplied = true;
    }
    return validateOptions(options, error);
}

// ==========================================
// Function: Resolve one dataset exposure list used by main-only mode
// Method: Prefer the single external/default list, otherwise derive the current
//         dataset's expo_<target>.list and normalize it to an absolute path.
// ==========================================
std::string resolveExposureList(const ProcessConfig::RuntimeOptions& options,
                                const ProcessConfig::DatasetSpec& dataset) {
    std::filesystem::path path;
    if (!options.expo_list.empty()) {
        path = options.expo_list;
    } else {
        if (options.output_root.empty()) {
            throw std::runtime_error(
                "--expo-list is absent and output-root cannot derive its default");
        }
        path = std::filesystem::path(options.output_root)
               / ("expo_" + dataset.target + ".list");
    }
    return std::filesystem::weakly_canonical(std::filesystem::absolute(path)).string();
}

// ==========================================
// Function: Format configured dataset defaults for help output
// Method: Join every target/prefix pair without mutating the configured list.
// ==========================================
std::string configuredDatasetsText() {
    if (ProcessConfig::DATASETS.empty()) {
        return "none";
    }
    std::string text;
    for (const ProcessConfig::DatasetSpec& dataset : ProcessConfig::DATASETS) {
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
    if (ProcessConfig::CONTAINS.empty()) {
        return "none (no token filter)";
    }
    std::string text;
    for (const std::string& token : ProcessConfig::CONTAINS) {
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
        << "  --run-init BOOL       Run initializer (default: "
        << (ProcessConfig::RUN_PROCESS_INIT ? "true" : "false") << ")\n"
        << "  --run-main BOOL       Run numerical pipeline (default: "
        << (ProcessConfig::RUN_PROCESS_MAIN ? "true" : "false") << ")\n"
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
        << ProcessConfig::EXISTING << ")\n"
        << "  --f77-max-path N      Generated path limit; zero disables (default: "
        << ProcessConfig::F77_MAX_PATH << ")\n"
        << "  --expo-list PATH      Single exposure list for main-only mode\n"
        << "  --help                Show this help\n"
        << "Options accept both --name value and --name=value. The first explicit "
           "--dataset or --contains replaces its configured list; repeats append.\n"
        << "Other duplicate scalar options use the last value. Main-only batches derive "
           "expo_<target>.list per dataset when --expo-list is omitted.\n"
        << "When both processes run, each process_init-generated absolute exposure-list "
           "path overrides external input for its dataset.\n";
}

}  // namespace

// ==========================================
// Function: Dispatch the integrated initializer and Fourier_Quad pipeline
// Method: Own MPI exactly once, parse one shared command line, process datasets
//         sequentially, and force chained processing to consume generated lists.
// ==========================================
int main(int argc, char* argv[]) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::my_id;
    int return_code = 0;

    ProcessConfig::RuntimeOptions options;
    std::string parse_error;
    const int local_parse_ok = parseCommandLine(argc, argv, options, parse_error) ? 1 : 0;
    int global_parse_ok = 0;
    MPI_Allreduce(&local_parse_ok, &global_parse_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (global_parse_ok == 0) {
        if (rank == 0) {
            std::cerr << "Argument error: "
                      << (parse_error.empty() ? "parsing failed on another rank" : parse_error)
                      << std::endl;
            printUsage(argv[0]);
        }
        return_code = 2;
    } else if (options.help_requested) {
        if (rank == 0) {
            printUsage(argv[0]);
        }
    } else {
        bool rng_initialized = false;
        for (std::size_t index = 0; index < options.datasets.size() && return_code == 0;
             ++index) {
            const ProcessConfig::DatasetSpec& dataset = options.datasets[index];
            if (rank == 0) {
                std::cout << "Dataset " << (index + 1) << "/" << options.datasets.size()
                          << ": target=" << dataset.target
                          << " prefix=" << dataset.prefix << std::endl;
            }

            std::string generated_exposure_list;
            if (options.run_process_init) {
                return_code = process_init(options, dataset, generated_exposure_list);
            }

            if (return_code == 0 && options.run_process_main) {
                std::string selected_exposure_list;
                if (options.run_process_init) {
                    selected_exposure_list = generated_exposure_list;
                    if (rank == 0) {
                        if (!options.expo_list.empty()) {
                            std::cout << "External exposure list overridden by process_init output: ";
                        } else {
                            std::cout << "process_main will use process_init output: ";
                        }
                        std::cout << selected_exposure_list << std::endl;
                    }
                    MPIScheduler::barrier();
                } else {
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
                    if (global_path_ok == 0) {
                        if (rank == 0) {
                            std::cerr << "Exposure-list argument error: "
                                      << (path_error.empty()
                                              ? "path resolution failed on another rank"
                                              : path_error)
                                      << std::endl;
                        }
                        return_code = 2;
                    }
                }

                if (return_code == 0 && !rng_initialized) {
                    const unsigned int rng_seed = NumericalRecipes::initializeRan1Seed(
                        rank, MPIScheduler::num_procs);
                    std::cout << "RNG_SEED rank seed: " << rank << " " << rng_seed << std::endl;
                    MPIScheduler::barrier();
                    rng_initialized = true;
                }
                if (return_code == 0) {
                    return_code = process_main(selected_exposure_list);
                }
            }

            if (return_code == 0 && index + 1 < options.datasets.size()) {
                MPIScheduler::barrier();
            }
        }
    }

    MPIScheduler::finalize();
    return return_code;
}
