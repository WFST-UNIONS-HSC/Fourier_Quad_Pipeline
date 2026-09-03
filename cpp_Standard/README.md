# cpp_Standard

Full C++17 Fourier_Quad pipeline with optional flat, mask, identity-astrometry,
external/hybrid PSF, and PCA/multi-scale branches.

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

Standard defaults to initialization, main, rearrangement, and FD enabled.
Review phase switches and paths before running. Fixed path defaults and output
layout names are centralized in `config/pathconfig.hpp`; CLI overrides workflow
paths. Scientific branches and thresholds in `config/LensingConfig.hpp` require
rebuilding.

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

Stage 9 keeps malformed paired shear rows fatal. On a float-parse failure it
now reports the failed zero/one-based column, token count and token, raw row
and hexadecimal tail, file/parser stream states, hostname, file stat metadata,
and a same-rank fresh-reopen comparison of the target row. The reopen is
diagnostic only: it never retries, repairs, or resumes the failed row, and the
pairing, quality cuts, calibration, and output schema are unchanged.

The focused `tests/CatalogCombinerLifecycleTest.cpp` regression covers the
diagnostic fields and representative finite float spellings. Build the main
program portably with `make -j4`; compile focused tests explicitly with the
same C++17 MPI wrapper and link settings. Locally, run the focused test binary
and `./Fourier_Quad_Pipe --help`. On a Linux cluster, load its MPI-enabled GCC,
CFITSIO, FFTW3, Eigen3, LAPACK, and BLAS modules, rebuild, and launch production
runs with the site's standard `srun` or `mpirun` command.

The default Stage-9 row now contains 18 external fields, `EXPO_NUM`, `ccD_NUM`,
and 25 pipeline fields (45 total). Regenerate Stage-9/rearr/FD products rather
than reading legacy 44-column catalogs with this build.
