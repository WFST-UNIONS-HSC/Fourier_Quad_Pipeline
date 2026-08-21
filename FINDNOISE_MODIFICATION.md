# Stage-3 local covariance noise-power update

This change permanently replaces the physical blank-noise-stamp path in C++ Standard and Lite.
Stage 3 now calls one `checkSourceAndEstimateNoisePower()` function for each retained source. The
function extracts one large local region, fits the existing source-local first-order plane once,
uses the same plane-subtracted residual for the recentered source and its surroundings, and
estimates a masked two-dimensional covariance from valid outer pixels on the source amplifier.

`NoiseCovariance` uses reusable FFTW plans and full `(2*N-1) x (2*N-1)` zero padding to compute
the residual autocorrelation numerator and mask pair counts. It retains signed short-lag
covariance, embeds the result in the existing `ns x ns` Fourier grid, and writes the real part of
the normalized transform without absolute values, clipping, square roots, random phases, or a
synthetic inverse transform.

## FITS data contract

- `*_source.fits`: real-space `ns x ns` source stamps, unchanged.
- `*_noise.fits`: Fourier-space, signed `ns x ns` local noise-power stamps.
- Stage 4 and Stage 6 copy the stored noise-power pixels directly into `noise_p`; only the source
  stamp is transformed. `processPowers()` continues the linear `source_p - noise_p` subtraction.

This is a producer/consumer contract change. Noise FITS files created by the old and new code
must not be mixed within one pipeline run.

## Compile-time controls

Both variants define the controls in `config/LensingConfig.hpp`:

| Parameter | Default | Meaning |
|---|---:|---|
| `noise_region_size` | 192 | Full local cutout side |
| `noise_inner_size` | 96 | Central square excluded from covariance |
| `noise_cov_max_lag` | 8 | Maximum retained signed x/y lag |
| `noise_cov_min_valid_pixels` | 4096 | Minimum valid outer pixels |
| `noise_cov_min_pair_fraction` | 0.50 | Minimum lag-pair count relative to zero lag |
| `noise_cov_sigma_ratio_min/max` | 0.80 / 1.25 | Allowed `sqrt(C(0,0))/sig` interval |
| `noise_cov_max_negative_fraction` | 0.25 | Catastrophic signed-negative-power QC limit |
| `noise_cov_imag_tolerance` | 1e-10 | Relative covariance-transform imaginary tolerance |

These values require representative-data runtime, rejection-rate, covariance-decay, and shear
bias validation before survey production tuning.

## Focused verification

`make test-noise-covariance` covers:

- FFT covariance versus direct pixel-pair covariance with no mask;
- a deterministic random mask;
- central exclusion, DQ holes, and same-amplifier clipping;
- anisotropic two-dimensional covariance;
- analytic covariance-to-power placement and normalization;
- preservation of signed negative modes;
- direct downstream copying of stored noise power without another FFT or magnitude square.

## Build environment

- Language standard: C++17.
- Local compiler/MPI: GCC `mpicxx` 15.2.0 and Open MPI 5.0.10.
- Local libraries: CFITSIO 4.6.4, FFTW 3.3.11, OpenBLAS BLAS/LAPACK 0.3.33, and Eigen 3.4.0.
- Cluster target: a Linux MPI toolchain with a C++17 compiler and compatible CFITSIO, FFTW,
  BLAS/LAPACK, and Eigen development modules. Use the site's supported module versions.

## Portable build and run

From either `cpp_Standard` or `cpp_Lite`:

```bash
make clean
make -j EIGEN_INCLUDE="$EIGEN_INCLUDE"
make test-noise-covariance EIGEN_INCLUDE="$EIGEN_INCLUDE"
./Fourier_Quad_Pipe --help
mpirun -np <ranks> ./Fourier_Quad_Pipe [options] [LEGACY_EXPO_LIST]
```

`EIGEN_INCLUDE` may be omitted when Eigen is in the compiler's default include path.

## Local WSL2 verification override

The 2026-08-21 verification used the existing Pixi compiler/library stack through
`STACK_PREFIX=/home/alatrion/.pixi/envs/base`, with `EIGEN_INCLUDE=/usr/include/eigen3` because
the Pixi base include directory did not contain Eigen. Standard and Lite each passed a clean full
build, `make test-noise-covariance`, all eight pre-existing Makefile regression targets, and
`./Fourier_Quad_Pipe --help`.
