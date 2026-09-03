# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for Gaia astrometry, no super-flat,
per-chip DQ masks, external sources, enabled deblending, local-polynomial PSF,
and no PCA. Alternate Standard branches are absent from this tree.

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

Pass `STACK_PREFIX` and, when needed, `EIGEN_INCLUDE` for nonstandard library
locations. The current Makefile exposes only `all` and `clean`; focused test
sources under `tests/` are compiled explicitly when needed.

Local focused verification uses the MPI C++ wrapper from GCC 15.2.0, with
CFITSIO 4.6.3 and FFTW3 3.3.10 available. The local full build uses Eigen3 from
`/usr/include/eigen3`; other sites must provide equivalent C++17 MPI, Eigen3,
LAPACK, and BLAS dependencies.

Lite defaults to initialization and main enabled, with rearrangement and FD
disabled. Its removed branches cannot be enabled by adding constants.
Fixed path defaults and output layout names are centralized in
`config/pathconfig.hpp`; runtime CLI values still override their compiled
workflow defaults.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. Its `--astrocat-output` directory
is independent of `LensingConfig::ASTROMETRY_CAT`; configure the consumer path
separately and set `LensingConfig::AstroCatType=2` before rebuilding when Stage
1 should read the generated tiles. `PathConfig::ASTROMETRY_TILE_PREFIX` defaults
to `astra_`; `PathConfig::SOURCE_CAT_TILE_PREFIX` defaults to `extern_` and is
shared by `process_extcat` and external-catalog lookup.

Stage 7 writes 24 fields and Stage 9 appends exposure chi-square. See the
[C++ guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).

The default Stage-9 row now contains 18 external fields, `EXPO_NUM`, `ccD_NUM`,
and 25 pipeline fields (45 total). Regenerate Stage-9/rearr/FD products rather
than reading legacy 44-column catalogs with this build.
