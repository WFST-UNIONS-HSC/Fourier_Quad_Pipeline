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
CFITSIO 4.6.4 and FFTW3 3.3.11 available. The local full build uses Eigen3 from
`/usr/include/eigen3`; other sites must provide equivalent C++17 MPI, Eigen3,
LAPACK, and BLAS dependencies.

Standard defaults to initialization, main, rearrangement, and FD enabled.
Review phase switches and paths before running. Fixed path defaults and output
layout names are centralized in `config/pathconfig.hpp`; CLI overrides workflow
paths. Scientific branches and thresholds in `config/LensingConfig.hpp` require
rebuilding.

Archive-format and detector naming conventions are compiled in
`config/InitConfig.hpp`: `ARCHIVE_SUFFIX` selects initializer inputs,
`CCDNUM_KEYWORD` names the DQ/main chip-number FITS keyword, and
`DQ_STEM_REPLACE_FROM`/`DQ_STEM_REPLACE_TO` map DQ archive stems to science
exposure stems. Their defaults remain `.fits.fz`, `CCDNUM`, and `ood` to `ooi`;
changing any of them requires `make clean` and a rebuild.

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

In the external-catalog path, Stage 9 counts every physical shear/orig line
with independent `getline` streams before opening the production readers. A
first mismatch triggers one fresh full recount; a persistent mismatch is
fatal, while a one-line shear catalog remains the legal header-only sentinel.
Matched files are then consumed for exactly the preflighted data-row count,
with both members of each pair read before row validation or science cuts.

Shear rows retain the fast stream-extraction parser and must provide all 24
floating fields. Stage 9 does not tokenize fields or add NaN/Inf
classification. The focused `tests/CatalogCombinerLifecycleTest.cpp`
regression covers both preflight attempts, fixed pairing, incomplete rows, and
representative finite float spellings. Build the main program portably with
`make -j4`; compile focused tests explicitly with the same C++17 MPI wrapper
and link settings. Locally, run the focused test binary and
`./Fourier_Quad_Pipe --help`. On a Linux cluster, load its MPI-enabled GCC,
CFITSIO, FFTW3, Eigen3, LAPACK, and BLAS modules, rebuild, and launch production
runs with the site's standard `srun` or `mpirun` command.

The default Stage-9 row now contains 18 external fields, `EXPO_NUM`, `ccD_NUM`,
and 25 pipeline fields (45 total). Regenerate Stage-9/rearr/FD products rather
than reading legacy 44-column catalogs with this build.
