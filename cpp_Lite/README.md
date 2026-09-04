# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for Gaia astrometry, no super-flat,
per-chip DQ masks, external sources, enabled deblending, local-polynomial PSF,
adaptive-pair PSF grouping, physical blank-noise stamps, and no PCA. Alternate
Standard and Lite selector branches are absent from this tree.

For ordinary use, edit `Initialize.hpp`, then build and run:

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 32 ./Fourier_Quad_Pipe
```

Runtime CLI options remain available for temporary or advanced overrides. For
example, `--source-cat PATH` overrides the main pipeline's external-catalog tile
input without changing `Initialize.hpp`.

Pass `STACK_PREFIX` and, when needed, `EIGEN_INCLUDE` for nonstandard library
locations. The current Makefile exposes only `all` and `clean`; focused test
sources under `tests/` are compiled explicitly when needed.

Local focused verification uses the MPI C++ wrapper from GCC 15.2.0, with
CFITSIO 4.6.3 and FFTW3 3.3.10 available. The local full build uses Eigen3 from
`/usr/include/eigen3`. A production cluster should load its site-provided C++17
MPI compiler plus CFITSIO, FFTW3, Eigen3, LAPACK, and BLAS modules; the portable
compile and run commands remain `make -j4` and
`mpirun -np <ranks> ./Fourier_Quad_Pipe`.

Lite defaults to initialization and main enabled, with rearrangement and FD
disabled. Its removed branches cannot be enabled by adding constants.
Ordinary instrument, workflow, dataset, path, catalog-schema, archive, and FITS
naming defaults are centralized in `Initialize.hpp`. `config/` remains the
internal compatibility and advanced-default layer; runtime CLI values still
override represented compiled workflow defaults.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. Its `--astrocat-output` directory
is independent of `Initialize::ASTROMETRY_CAT`, and Lite Stage 1 always reads
the configured repartitioned one-degree catalog layout. Likewise,
`--extcat-output` controls only the `process_extcat` producer, while
`Initialize::SOURCE_CAT_DEFAULT` or `--source-cat` controls the `process_main`
consumer. Producer output defaults intentionally remain empty in the advanced
config layer, so enabling either producer requires its output CLI option.
`Initialize::ASTROMETRY_TILE_PREFIX` defaults to `astra_`, and
`Initialize::SOURCE_CAT_TILE_PREFIX` defaults to `extern_`.

Stage 7 writes 24 fields and Stage 9 appends exposure chi-square. See the
[C++ guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).

Stage 9 counts every physical shear/orig line with independent `getline`
streams before opening the production readers. A first mismatch triggers one
fresh full recount; a persistent mismatch is fatal, while a one-line shear
catalog remains the legal header-only sentinel. Matched files are consumed for
exactly the preflighted data-row count, with both members of each pair read
before row validation or science cuts. Shear rows retain the fast
stream-extraction parser, must provide all 24 floating fields, and are not
tokenized for separate NaN/Inf classification. The focused
`tests/CatalogCombinerLifecycleTest.cpp` regression covers both preflight
attempts, fixed pairing, and incomplete rows.

The default Stage-9 row now contains 18 external fields, `EXPO_NUM`, `ccD_NUM`,
and 25 pipeline fields (45 total). Regenerate Stage-9/rearr/FD products rather
than reading legacy 44-column catalogs with this build.
