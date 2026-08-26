# cpp_Lite

The frozen-branch simplified C++17 Fourier_Quad pipeline build (`PSFRecons`
removed). See `REFACTOR_NOTES.md` for the Lite change log.

The complete C++ pipeline guide lives in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference is
[`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Stage-3 outer-noise plane fitting

For covariance noise products (`NstampType=2`), Stage 3 fits the source and
noise residual plane from the configurable square shell between
`noise_region_size` and `noise_inner_size`. It rejects masked, non-finite,
out-of-chip, and opposite-amplifier samples. This retained Stage-3 path does not
restore any Lite-deleted PSF/PCA branch. The synthetic regression is kept
independent of the Makefile and can be compiled with the production plane
solver as follows:

```bash
mpicxx -O2 -std=c++17 -Wall -Wextra -ffunction-sections -fdata-sections \
  -Iinclude -Iconfig -Iinclude/process_main -Isrc/process_main \
  -I"${STACK_PREFIX}/include" -I"${EIGEN_INCLUDE}" \
  tests/NoisePlaneFitTest.cpp src/process_main/UniversalUtils.cpp \
  src/process_main/LinearSolve.cpp -Wl,--gc-sections \
  -L"${STACK_PREFIX}/lib" -Wl,-rpath,"${STACK_PREFIX}/lib" \
  -llapack -lblas -lm -o /tmp/NoisePlaneFitTest
/tmp/NoisePlaneFitTest
```

## Stage-5 PSF star-selection redesign

Stage 5 now applies a positive signed central-chi-window quality gate, reads the
matched image positions already stored in each chip's Gaia `_astro.dat`,
estimates a Gaia-assisted exposure-wide FWHM stellar locus, and performs only
same-chip Fourier comparisons. `LensingConfig::PsfGroupingType` selects either
the legacy chi-threshold graph (`1`) or exact mutual-KNN graph (`2`); both paths
share the same exposure-pooled per-star `minChi` cut and
main/eligible-secondary component selection. Type 2 rebuilds its exact top-K
lists only among candidates that survive that shared `minChi` cut, so rejected
neighbours cannot occupy stale KNN slots; Type 1 retains its existing private
threshold sample and graph behavior. Secondary components must pass both the
configured relative-size and Gaia-count conditions. The former candidate-count-
squared chi matrix is gone.

The retained stars receive one analytic leave-one-out PRESS pass. Chips reuse
their initial normalized local-polynomial fit when no star is removed and refit
exactly once when the retained set changes. Lite's `msshape_*` diagnostics are
produced from the final analytic LOO model. This port deliberately does not add
the Standard-only trivial-astrometry, PCA/multi-scale, residual-reconstruction,
or hybrid PSF branches that Lite removed.

All new scientific constants are compile-time values in
`config/LensingConfig.hpp`; changing them requires rebuilding. The focused
selection tests are:

```bash
make test-psf-star-selection CXX=mpicxx \
     STACK_PREFIX=/path/to/dependency-prefix \
     EIGEN_INCLUDE=/path/to/eigen3
```

Synthetic tests cover the shared chi window and quality gate, FWHM/Gaia peak
selection, `_astro.dat` parsing and nearest matching, streaming top-K and mutual
components, survivor-only KNN slot refill after `minChi`, the shared secondary-
group policy, non-square state storage, and analytic LOO equivalence to explicit
leave-one-out refits. Representative real exposures are still required to
inspect PRESS versus brightness/SNR/FWHM and to benchmark Stage-5 wall time.

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
FITS open, subtracts the Stage-1 background model after the retained Lite
preprocessing state, and reconstructs the sigma map from the header. Missing or
malformed metadata is a chip failure; there is no legacy pixel-metadata
fallback. Lite's frozen branches remain unchanged.

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
