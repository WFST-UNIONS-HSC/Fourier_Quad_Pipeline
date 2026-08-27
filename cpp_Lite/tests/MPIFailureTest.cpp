#include "general/MPIScheduler.hpp"
#include "process_main/MPIFailure.hpp"

#include <cstdlib>

// ==========================================
// Function: Exercise MPI-wide fatal termination from one worker
// Method: Abort rank one while rank zero waits in a collective; success is a fast nonzero job exit.
// ==========================================
int main(int argc, char** argv) {
    MPIScheduler::init(argc, argv);
    const int rank = MPIScheduler::state.rank;
    if (rank == 1) {
        MPIFailure::abortWorld(
            "intentional MPI failure test", "rank-one synthetic missing input");
    }
    MPIScheduler::barrier();
    MPIScheduler::finalize();
    return EXIT_SUCCESS;
}
