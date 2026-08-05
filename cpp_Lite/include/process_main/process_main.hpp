#ifndef PROCESS_MAIN_PROCESS_MAIN_HPP
#define PROCESS_MAIN_PROCESS_MAIN_HPP

#include <string>

// ==========================================
// Function: Run the numerical Fourier_Quad pipeline on one exposure list
// Method: Load and broadcast the list, then execute the configured MPI stages
//         without owning MPI initialization or finalization.
// ==========================================
int process_main(const std::string& exposure_list);

#endif  // PROCESS_MAIN_PROCESS_MAIN_HPP
