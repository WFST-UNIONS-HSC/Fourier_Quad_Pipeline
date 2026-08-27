#include "process_main/MPIFailure.hpp"

#include <mpi.h>

#include <cstdlib>

// ==========================================
// Function: Exercise MPI-wide fatal termination from one worker
// Method: Abort rank one while rank zero waits in a collective; success is a fast nonzero job exit.
// ==========================================
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 1) {
        MPIFailure::abortWorld(
            "intentional MPI failure test", "rank-one synthetic missing input");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
