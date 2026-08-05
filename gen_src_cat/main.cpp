#include "process_extcat/process_extcat.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ==========================================
// Structure: Track repeatable CLI override state
// Method: Clear configured filename tokens only when the first explicit filter appears.
// ==========================================
struct ParserState {
    bool contains_option_seen = false;
};

// ==========================================
// Function: Parse one command-line boolean
// Method: Accept the same textual and numeric forms used by the integrated pipeline CLI.
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
// Function: Parse one positive unsigned integer
// Method: Require full decimal consumption and reject zero or overflow.
// ==========================================
bool parsePositiveInteger(const std::string& value, std::uint64_t& parsed) {
    if (value.empty()
        || !std::all_of(value.begin(), value.end(), [](unsigned char character) {
               return character >= '0' && character <= '9';
           })) {
        return false;
    }
    std::size_t consumed = 0;
    unsigned long long number = 0;
    try {
        number = std::stoull(value, &consumed);
    } catch (const std::exception&) {
        return false;
    }
    if (consumed != value.size() || number == 0) {
        return false;
    }
    parsed = static_cast<std::uint64_t>(number);
    return true;
}

// ==========================================
// Function: Parse the input delimiter name
// Method: Map one stable lowercase CLI value to the public configuration enum.
// ==========================================
bool parseDelimiter(const std::string& value, ProcessExtcat::Delimiter& delimiter) {
    if (value == "auto") {
        delimiter = ProcessExtcat::Delimiter::Auto;
    } else if (value == "whitespace") {
        delimiter = ProcessExtcat::Delimiter::Whitespace;
    } else if (value == "comma") {
        delimiter = ProcessExtcat::Delimiter::Comma;
    } else if (value == "tab") {
        delimiter = ProcessExtcat::Delimiter::Tab;
    } else {
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the input header policy
// Method: Map auto, present, or absent to the public configuration enum.
// ==========================================
bool parseHeaderMode(const std::string& value, ProcessExtcat::HeaderMode& mode) {
    if (value == "auto") {
        mode = ProcessExtcat::HeaderMode::Auto;
    } else if (value == "present") {
        mode = ProcessExtcat::HeaderMode::Present;
    } else if (value == "absent") {
        mode = ProcessExtcat::HeaderMode::Absent;
    } else {
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the malformed-row policy
// Method: Accept only fail-fast or explicit skip behavior.
// ==========================================
bool parseMalformedPolicy(const std::string& value,
                          ProcessExtcat::MalformedPolicy& policy) {
    if (value == "fail") {
        policy = ProcessExtcat::MalformedPolicy::Fail;
    } else if (value == "skip") {
        policy = ProcessExtcat::MalformedPolicy::Skip;
    } else {
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the existing-output policy
// Method: Accept only fail-fast or transactional overwrite behavior.
// ==========================================
bool parseExistingPolicy(const std::string& value,
                         ProcessExtcat::ExistingPolicy& policy) {
    if (value == "fail") {
        policy = ProcessExtcat::ExistingPolicy::Fail;
    } else if (value == "overwrite") {
        policy = ProcessExtcat::ExistingPolicy::Overwrite;
    } else {
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse an explicit canonical column projection
// Method: Convert exactly 18 comma-separated one-based indices to zero-based indices.
// ==========================================
bool parseColumns(const std::string& value,
                  std::array<std::size_t, ProcessExtcat::kCanonicalColumnCount>& columns) {
    std::istringstream input(value);
    std::string token;
    std::vector<std::size_t> parsed;
    while (std::getline(input, token, ',')) {
        std::uint64_t index = 0;
        if (!parsePositiveInteger(token, index)
            || index - 1 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            return false;
        }
        parsed.push_back(static_cast<std::size_t>(index - 1));
    }
    if (parsed.size() != ProcessExtcat::kCanonicalColumnCount) {
        return false;
    }
    for (std::size_t index = 0; index < parsed.size(); ++index) {
        columns[index] = parsed[index];
    }
    return true;
}

// ==========================================
// Function: Read the value attached to one named option
// Method: Support both --name=value and --name value without accepting missing values.
// ==========================================
bool optionValue(const std::string& argument,
                 int& index,
                 int argc,
                 char** argv,
                 std::string& name,
                 std::string& value,
                 std::string& error) {
    const std::size_t equals = argument.find('=');
    if (equals != std::string::npos) {
        name = argument.substr(0, equals);
        value = argument.substr(equals + 1);
        if (value.empty()) {
            error = "missing value after " + name;
            return false;
        }
        return true;
    }
    name = argument;
    if (index + 1 >= argc) {
        error = "missing value after " + name;
        return false;
    }
    ++index;
    value = argv[index];
    return true;
}

// ==========================================
// Function: Apply one named command-line option
// Method: Update the standalone config while preserving repeatable OR-matched filters.
// ==========================================
bool applyOption(const std::string& name,
                 const std::string& value,
                 ProcessExtcat::Config& config,
                 ParserState& state,
                 std::string& error) {
    if (name == "--input-dir") {
        config.input_directory = value;
    } else if (name == "--output-dir") {
        config.output_directory = value;
    } else if (name == "--contains") {
        if (value.empty()) {
            error = "--contains must not be empty";
            return false;
        }
        if (!state.contains_option_seen) {
            config.filename_tokens.clear();
            state.contains_option_seen = true;
        }
        config.filename_tokens.push_back(value);
    } else if (name == "--recursive") {
        if (!parseBoolean(value, config.recursive)) {
            error = "--recursive must be true, false, 1, 0, on, or off";
            return false;
        }
    } else if (name == "--delimiter") {
        if (!parseDelimiter(value, config.delimiter)) {
            error = "--delimiter must be auto, whitespace, comma, or tab";
            return false;
        }
    } else if (name == "--header") {
        if (!parseHeaderMode(value, config.header_mode)) {
            error = "--header must be auto, present, or absent";
            return false;
        }
    } else if (name == "--columns") {
        if (!parseColumns(value, config.input_columns)) {
            error = "--columns must contain exactly 18 positive one-based indices";
            return false;
        }
        config.use_explicit_columns = true;
    } else if (name == "--chunk-mib") {
        std::uint64_t mebibytes = 0;
        if (!parsePositiveInteger(value, mebibytes)
            || mebibytes > std::numeric_limits<std::uint64_t>::max() / (1024ULL * 1024ULL)) {
            error = "--chunk-mib must be a positive integer within uint64 range";
            return false;
        }
        config.chunk_bytes = mebibytes * 1024ULL * 1024ULL;
    } else if (name == "--malformed") {
        if (!parseMalformedPolicy(value, config.malformed_policy)) {
            error = "--malformed must be fail or skip";
            return false;
        }
    } else if (name == "--existing") {
        if (!parseExistingPolicy(value, config.existing_policy)) {
            error = "--existing must be fail or overwrite";
            return false;
        }
    } else {
        error = "unknown option: " + name;
        return false;
    }
    return true;
}

// ==========================================
// Function: Parse the standalone process_extcat command line
// Method: Accept named options in any order and reject positional arguments.
// ==========================================
bool parseCommandLine(int argc,
                      char** argv,
                      ProcessExtcat::Config& config,
                      bool& help_requested,
                      std::string& error) {
    ParserState state;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            help_requested = true;
            continue;
        }
        if (argument.rfind("--", 0) != 0) {
            error = "unexpected positional argument: " + argument;
            return false;
        }
        std::string name;
        std::string value;
        if (!optionValue(argument, index, argc, argv, name, value, error)
            || !applyOption(name, value, config, state, error)) {
            return false;
        }
    }
    return true;
}

// ==========================================
// Function: Print standalone process_extcat usage
// Method: Document the portable CLI and the fixed canonical output contract.
// ==========================================
void printUsage(const char* executable) {
    std::cout
        << "Usage: mpirun -np <N> " << executable << " --input-dir PATH --output-dir PATH [options]\n"
        << "\n"
        << "Options:\n"
        << "  --contains TEXT       Repeatable basename substring, matched with OR; default: all files\n"
        << "  --recursive BOOL      Recurse below input directory (default: true)\n"
        << "  --delimiter MODE      auto, whitespace, comma, or tab (default: auto)\n"
        << "  --header MODE         auto, present, or absent (default: auto)\n"
        << "  --columns LIST        18 comma-separated one-based input indices in canonical output order\n"
        << "  --chunk-mib N         Maximum nominal byte-range task size (default: 64)\n"
        << "  --malformed POLICY    fail or skip (default: fail)\n"
        << "  --existing POLICY     fail or overwrite generated tiles (default: fail)\n"
        << "  --help, -h            Show this help text\n"
        << "\n"
        << "Output is always one-degree des_y6_RA_*_Dec_*.dat tiles with the canonical\n"
        << "18-column commented-header schema consumed by the Fourier_Quad pipeline.\n";
}

}  // namespace

// ==========================================
// Function: Run the standalone MPI external-catalog tiler
// Method: Own MPI only in this wrapper so process_extcat remains directly embeddable.
// ==========================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ProcessExtcat::Config config;
    bool help_requested = false;
    std::string parse_error;
    const int local_parse_ok = parseCommandLine(
        argc, argv, config, help_requested, parse_error) ? 1 : 0;
    int global_parse_ok = 0;
    MPI_Allreduce(&local_parse_ok, &global_parse_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    int return_code = 0;
    if (global_parse_ok == 0) {
        if (rank == 0) {
            std::cerr << "Argument error: "
                      << (parse_error.empty() ? "parsing failed on another rank" : parse_error)
                      << "\n";
            printUsage(argv[0]);
        }
        return_code = 2;
    } else if (help_requested) {
        if (rank == 0) {
            printUsage(argv[0]);
        }
    } else {
        return_code = process_extcat(config, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return return_code;
}
