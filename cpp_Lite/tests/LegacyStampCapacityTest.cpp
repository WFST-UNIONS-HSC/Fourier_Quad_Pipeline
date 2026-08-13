#include "FitsIO.hpp"

#include <fitsio.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

namespace {

// ==========================================
// Function: Stop the legacy-stamp test when one invariant fails
// Method: Print a focused diagnostic and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "LegacyStampCapacity test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Require a FITS product to retain two image axes
// Method: Read NAXIS through CFITSIO without relying on the production reader.
// ==========================================
void requireNaxis2(const std::string& filename) {
    fitsfile* file = nullptr;
    int status = 0;
    int naxis = 0;
    fits_open_file(&file, filename.c_str(), READONLY, &status);
    fits_get_img_dim(file, &naxis, &status);
    fits_close_file(file, &status);
    require(status == 0 && naxis == 2, "stamp product must remain NAXIS=2");
}

// ==========================================
// Function: Verify one above/below-threshold legacy mosaic round trip
// Method: Use tall 1x64 stamps to cross the obsolete 7000-row guard cheaply.
// ==========================================
void testCount(const std::string& filename, int count) {
    constexpr int nsx = 1;
    constexpr int nsy = 64;
    constexpr int row_length = 40;
    const int n1 = nsx * row_length;
    const int n2 = nsy * (count / row_length + 1);

    std::vector<float> input(static_cast<std::size_t>(count) * nsx * nsy);
    for (std::size_t index = 0; index < input.size(); ++index) {
        input[index] = static_cast<float>(static_cast<int>(index % 251) - 125);
    }

    require(FitsIO::writeStamps(
                count, 1, count, nsx, nsy, input, n1, n2, filename),
            "legacy stamp write must succeed");
    requireNaxis2(filename);

    std::vector<float> output;
    require(FitsIO::readStamps(
                count, 1, count, nsx, nsy, output, n1, n2, filename),
            "legacy stamp read must succeed");
    require(output == input, "legacy NAXIS=2 stamp round trip must be exact");
}

// ==========================================
// Function: Verify selected-stamp packing beyond the former exposure limit
// Method: Select every live stamp, write through writeStamps2, and round-trip
//         the unchanged two-dimensional mosaic.
// ==========================================
void testSelectedCount(const std::string& filename, int count) {
    constexpr int nsx = 1;
    constexpr int nsy = 64;
    constexpr int row_length = 40;
    const int n1 = nsx * row_length;
    const int n2 = nsy * (count / row_length + 1);

    std::vector<float> input(static_cast<std::size_t>(count) * nsx * nsy);
    for (std::size_t index = 0; index < input.size(); ++index) {
        input[index] = static_cast<float>(static_cast<int>(index % 197) - 98);
    }
    const std::vector<int> selected(static_cast<std::size_t>(count), 1);

    require(FitsIO::writeStamps2(
                count, count, nsx, nsy, input, selected, 1, n1, n2,
                filename),
            "selected legacy stamp write must succeed");
    requireNaxis2(filename);

    std::vector<float> output;
    require(FitsIO::readStamps(
                count, 1, count, nsx, nsy, output, n1, n2, filename),
            "selected legacy stamp read must succeed");
    require(output == input,
            "selected legacy NAXIS=2 stamp round trip must be exact");
}

}  // namespace

// ==========================================
// Function: Run the legacy NAXIS=2 capacity regression suite
// Method: Exercise boundary-adjacent and large live stamp counts in one temporary file.
// ==========================================
int main() {
    const std::filesystem::path filename =
        std::filesystem::temp_directory_path()
        / ("fq_legacy_stamp_" + std::to_string(::getpid()) + ".fits");
    const int counts[] = {3999, 4000, 4001, 5000, 8000};
    for (int count : counts) {
        testCount(filename.string(), count);
    }
    const int selected_counts[] = {4999, 5000, 5001, 7000};
    for (int count : selected_counts) {
        testSelectedCount(filename.string(), count);
    }
    std::error_code error;
    std::filesystem::remove(filename, error);
    std::cout << "LegacyStampCapacity tests passed\n";
    return EXIT_SUCCESS;
}
