# cpp_Lite

The frozen-branch simplified C++17 Fourier_Quad pipeline build (`PSFRecons`
removed). See `REFACTOR_NOTES.md` for the Lite change log.

The complete C++ pipeline guide lives in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference is
[`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Validated build and regression checks

The v1.3.1 capacity and failure-handling update was validated under WSL2 with
C++17 GCC 15.2.0, Open MPI 5.0.10, Eigen3, CFITSIO, FFTW3, LAPACK, and BLAS.
Build with the MPI compiler wrapper and the dependency prefixes available on
the target system:

```bash
make CXX=mpicxx STACK_PREFIX=/path/to/dependency-prefix EIGEN_INCLUDE=/path/to/eigen3
```

The focused regression targets cover Stage-1 norm gating, Stage-9 catalog row
counts, dynamic PSF state, FITS stamp counts beyond the legacy limits, and
dynamic astrometry catalog sizes:

```bash
make test-universalblock test-catalog-row-count test-psf-model-state \
     test-legacy-stamp-capacity test-astrometry-dynamic-capacity \
     test-catalog-lifecycle test-exposure-runtime-sizing \
     test-mpi-failure \
     CXX=mpicxx STACK_PREFIX=/path/to/dependency-prefix \
     EIGEN_INCLUDE=/path/to/eigen3 MPIRUN=mpirun
```

The MPI integration target treats a quick nonzero two-rank termination as
success after one worker reports the injected fatal error; a zero status or
ten-second timeout fails the target.
