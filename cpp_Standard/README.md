# cpp_Standard

The full C++17 Fourier_Quad pipeline build (includes PCA `PSFRecons`).

The complete C++ pipeline guide - source structure, pipeline stages,
configuration, building, run modes, initializer output layout, Docker, and HPC -
now lives in [`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference
is [`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

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

## Stage-1 norm FITS coefficient contract

Stage 1 writes its final background polynomial coefficients and sigma-plane
coefficients into the `norm.fits` primary header. `CCD_split=1` uses the
long-string keywords `BGCO` and `SIGCO`; `CCD_split=2` uses `BG1CO`, `BG2CO`,
`SIG1CO`, and `SIG2CO`. Stage 3 reads the image and these keywords through one
FITS open, subtracts the Stage-1 background model after optional flat
correction, and reconstructs the sigma map from the header. Missing or
malformed metadata is a chip failure; there is no legacy pixel-metadata
fallback.

Local WSL2 verification used GCC 15.2.0, Open MPI 5.0.10, CFITSIO, FFTW3,
LAPACK, BLAS, and Eigen3. The portable build command remains:

```bash
make CXX=mpicxx STACK_PREFIX=/path/to/dependency-prefix \
     EIGEN_INCLUDE=/path/to/eigen3
./Fourier_Quad_Pipe --help
```

The local override used for this checkout was
`CXX=/home/alatrion/.pixi/bin/mpicxx`,
`STACK_PREFIX=/home/alatrion/.pixi/envs/base`, and
`EIGEN_INCLUDE=/usr/include/eigen3`. The production target should provide
equivalent compiler/library modules on Linux HPC.
