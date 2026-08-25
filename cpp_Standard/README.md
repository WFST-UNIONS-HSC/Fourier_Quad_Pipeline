# cpp_Standard

The full C++17 Fourier_Quad pipeline build (includes PCA `PSFRecons`).

The complete C++ pipeline guide - source structure, pipeline stages,
configuration, building, run modes, initializer output layout, Docker, and HPC -
now lives in [`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference
is [`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Stage-5 PSF star-selection redesign

Stage 5 now applies a positive signed central-chi-window quality gate, reads the
matched image positions already stored in each chip's `_astro.dat`, estimates a
Gaia-assisted exposure-wide FWHM stellar locus, and performs only same-chip
Fourier comparisons. `LensingConfig::PsfGroupingType` selects either the legacy
chi-threshold graph (`1`) or exact mutual-KNN graph (`2`); both paths share the
same exposure-pooled per-star `minChi` cut and main/eligible-secondary component
selection. Type 2 rebuilds its exact top-K lists only among candidates that
survive that shared `minChi` cut, so rejected neighbours cannot occupy stale KNN
slots; Type 1 retains its existing private threshold sample and graph behavior.
Secondary components must pass both the configured relative-size and Gaia-count
conditions. The former candidate-count-squared chi matrix is gone.

The retained stars receive one analytic leave-one-out PRESS pass. Chips reuse
their initial polynomial fit when no star is removed and refit exactly once when
the retained set changes. For the local-polynomial branch, `msshape_*` and the
optional PCA residual stamps are produced from the final analytic LOO model and
residual. The optional hybrid branch shares the new selection/PRESS and cached
native-coordinate polynomial fit, while retaining its separate very-local model
diagnostic because the polynomial hat-diagonal formula alone does not define an
exact leave-one-out hybrid interpolation.

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
FITS open, subtracts the Stage-1 background model after optional flat
correction, and reconstructs the sigma map from the header. Missing or
malformed metadata is a chip failure; there is no legacy pixel-metadata
fallback.

When `include_FLAT=1`, Standard Stage 3 reuses the Stage-1 flat filename
`<flat_path>/flat_<two-digit-chip>_weight.fits` (the compiled
`LensingConfig::FLAT_PATH` in this repository). Flat read failures and dimension
mismatches are fatal chip errors. For each pixel, `flat < 0.5` masks the weight;
otherwise the science array is multiplied by the flat before the recorded
background model is subtracted.

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
