# Stage-3 local covariance noise-power update

This change permanently replaces the physical blank-noise-stamp path in C++ Standard and Lite.
Stage 3 now calls one `checkSourceAndEstimateNoisePower()` function for each retained source. The
function extracts one large local region, fits the existing source-local first-order plane once,
uses the same plane-subtracted residual for the recentered source and its surroundings, and
estimates a masked two-dimensional covariance from valid outer pixels on the source amplifier.

`NoiseCovariance` uses reusable FFTW plans and explicit factor-2 zero padding to compute the
residual autocorrelation numerator and mask pair counts. For the 192-pixel production region this
is a `384 x 384` FFT, with a compile-time check that the side remains at least `2*N-1`; workspace
reuse is keyed by both region and padded sides. It retains signed short-lag covariance, embeds the
result in the existing `ns x ns` Fourier grid, and writes the normalized real transform without
absolute values, clipping, square roots, random phases, or a synthetic inverse transform.

After recentering, a source is rejected when its final `ns x ns` stamp straddles the midpoint of a
two-amplifier chip. Stamps ending at or starting on the boundary remain valid.

For an accepted source whose final stamp still contains peripheral masked pixels, Stage 3 now waits
until covariance and signed-power QC pass, clips negative modes only in a temporary synthesis PSD,
and renormalizes that PSD so its sum equals the measured `C(0,0)`. An ifftshifted Hermitian Gaussian
field is inverse transformed and only same-coordinate `weight == 0` pixels are replaced. Clean stamps
return before FFT allocation or RNG use; valid pixels and stored signed noise power remain unchanged.
Decoration failure rejects the source and never falls back to independent white noise.

## FITS data contract

- `*_source.fits`: real-space `ns x ns` source stamps; retained peripheral masked pixels contain a
  local correlated-noise realization instead of independent white Gaussian draws.
- `*_noise.fits`: Fourier-space, signed `ns x ns` local noise-power stamps.
- Stage 4 and Stage 6 copy stored noise power directly. Source stamps always use raw
  `getPower(..., 0)`, including the source-only SNR_F diagnostic path.
- Measurement power follows one shared order: subtract stored noise, apply `star_smooth` or
  `gal_smooth` to corrected power, then subtract its outer-edge mean. Noise power is never smoothed
  separately. Stars are regularized after this sequence.
- Log smoothing is signed-safe: it retains the legacy span offset when sufficient and otherwise
  raises the offset above `-Pmin` by a scale-dependent epsilon before log/smooth/exp inversion.

This is a producer/consumer contract change. Noise FITS files created by the old and new code
must not be mixed within one pipeline run.

## Compile-time controls

Both variants define the controls in `config/LensingConfig.hpp`:

| Parameter | Default | Meaning |
|---|---:|---|
| `noise_region_size` | 192 | Full local cutout side |
| `noise_inner_size` | 96 | Central square excluded from covariance |
| `noise_cov_padding_factor` | 2.0 | Linear-autocorrelation FFT padding factor |
| `noise_cov_fft_size` | 384 | Derived padded FFT side |
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
- the production 192-to-384 padding derivation, undersized-padding rejection, and cache re-keying;
- source stamps immediately beside or crossing the amplifier boundary.

`make test-power-processing` covers:

- finite signed-log smoothing for positive, zero, and negative corrected modes;
- positive-only compatibility with the former logarithmic smoother;
- raw source FFT followed by subtract, configured smooth mode 0/1/2, and edge subtraction;
- rejection of the former smoothing-before-subtraction order;
- zero post-subtraction outer-edge mean and final star regularization.

`make test-correlated-decoration` covers:

- bitwise no-mask and unmasked-pixel preservation, including no RNG consumption on the fast path;
- finite same-coordinate masked replacement and hard failure without white-noise fallback;
- source-core, connected-neighbor, and peripheral-mask protection in `markSource`;
- signed stored-power immutability and clipped synthesis-PSD renormalization;
- ensemble Fourier-power, real-space variance, and anisotropy closure.

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
make test-power-processing EIGEN_INCLUDE="$EIGEN_INCLUDE"
make test-correlated-decoration EIGEN_INCLUDE="$EIGEN_INCLUDE"
./Fourier_Quad_Pipe --help
mpirun -np <ranks> ./Fourier_Quad_Pipe [options] [LEGACY_EXPO_LIST]
```

`EIGEN_INCLUDE` may be omitted when Eigen is in the compiler's default include path.

## Local WSL2 verification override

The 2026-08-21 verification used the existing Pixi compiler/library stack through
`STACK_PREFIX=/home/alatrion/.pixi/envs/base`, with `EIGEN_INCLUDE=/usr/include/eigen3` because
the Pixi base include directory did not contain Eigen. Standard and Lite each passed a clean full
build, `make test-noise-covariance`, `make test-power-processing`, all eight pre-existing Makefile
regression targets, `make test-correlated-decoration`, and `./Fourier_Quad_Pipe --help`.
