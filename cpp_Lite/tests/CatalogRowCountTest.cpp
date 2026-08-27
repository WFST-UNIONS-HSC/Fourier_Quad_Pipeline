#include "process_main/CatalogRowCount.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>

#include <sys/wait.h>
#include <unistd.h>

namespace {

// ==========================================
// Class: Own one temporary physical-row catalog
// Method: Write an optional header plus an exact data-line count and remove it afterward.
// ==========================================
class TemporaryCatalog {
public:
    // ==========================================
    // Function: Create one temporary physical-row catalog
    // Method: Write the optional header and requested exact data-line count.
    // ==========================================
    TemporaryCatalog(const std::string& label, std::size_t data_rows,
                     bool write_header = true)
        : path_(std::filesystem::temp_directory_path()
                / ("fq_catalog_rows_" + label + "_"
                   + std::to_string(::getpid()) + ".dat")) {
        std::ofstream output(path_);
        if (write_header) {
            output << "header\n";
        }
        for (std::size_t row = 0; row < data_rows; ++row) {
            output << row << '\n';
        }
    }

    // ==========================================
    // Function: Remove one temporary catalog
    // Method: Clean the file without propagating teardown errors.
    // ==========================================
    ~TemporaryCatalog() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    // ==========================================
    // Function: Return the temporary catalog path
    // Method: Serialize the stored filesystem path for production helpers.
    // ==========================================
    std::string path() const { return path_.string(); }

private:
    std::filesystem::path path_;
};

// ==========================================
// Function: Stop the row-count test when one invariant fails
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "CatalogRowCount test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Require one fail-fast catalog operation to terminate a child
// Method: Fork an isolated caller and accept only a nonzero exit or signal.
// ==========================================
template <typename Operation>
void requireAbort(Operation operation, const std::string& message) {
    const pid_t child = ::fork();
    require(child >= 0, "fork failed for fatal case");
    if (child == 0) {
        operation();
        std::_Exit(EXIT_SUCCESS);
    }

    int status = 0;
    require(::waitpid(child, &status, 0) == child,
            "waitpid failed for fatal case");
    require((WIFEXITED(status) && WEXITSTATUS(status) != EXIT_SUCCESS)
                || WIFSIGNALED(status),
            message);
}

// ==========================================
// Function: Exercise equal, asymmetric, zero-row, and missing-header catalog cases
// Method: Use production helpers for every physical-row decision.
// ==========================================
void testRowCountCases() {
    const TemporaryCatalog rows100("rows100", 100);
    const TemporaryCatalog rows99("rows99", 99);
    const TemporaryCatalog zero_a("zero_a", 0);
    const TemporaryCatalog zero_b("zero_b", 0);
    const TemporaryCatalog no_header("no_header", 0, false);

    require(CatalogCombiner::Internal::countCatalogDataRows(rows100.path()) == 100,
            "catalog must contain 100 data rows");
    CatalogCombiner::Internal::requireMatchingCatalogDataRows(
        rows100.path(), rows100.path());
    CatalogCombiner::Internal::requireMatchingCatalogDataRows(
        zero_a.path(), zero_b.path());

    requireAbort(
        [&]() {
            CatalogCombiner::Internal::requireMatchingCatalogDataRows(
                rows100.path(), rows99.path());
        },
        "100/99 rows must abort");
    requireAbort(
        [&]() {
            CatalogCombiner::Internal::requireMatchingCatalogDataRows(
                rows99.path(), rows100.path());
        },
        "99/100 rows must abort");
    requireAbort(
        [&]() {
            CatalogCombiner::Internal::countCatalogDataRows(no_header.path());
        },
        "missing header must abort");
}

}  // namespace

// ==========================================
// Function: Run focused Stage-9 physical-row regression cases
// Method: Execute all synthetic catalog pairs and report one success line.
// ==========================================
int main() {
    testRowCountCases();
    std::cout << "CatalogRowCount tests passed\n";
    return EXIT_SUCCESS;
}
