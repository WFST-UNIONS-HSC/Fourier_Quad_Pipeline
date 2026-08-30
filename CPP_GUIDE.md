# Fourier_Quad C++ Pipeline Guide

This guide describes the C++ programs in this repository. They use CLI
overrides and compiled defaults; unlike the separate `Fourier_Quad_Cpp`
repository, this version has no INI configuration layer.

> 中文版：[CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)

## Variants and layout

[`cpp_Standard`](cpp_Standard/) (C++ Standard) retains optional flat, mask,
identity astrometry, external/hybrid PSF, and PCA branches.

[`cpp_Lite`](cpp_Lite/) (C++ Lite) physically removes those alternatives and
keeps Gaia astrometry, per-chip DQ
masks, External source catalog matching, deblending, local-polynomial PSF, and
no PCA. Because C++ Lite fixes the per-chip DQ branch, its runs require DQ masks;
C++ Standard can select a configuration that does not read them.

Each variant contains `main.cpp`, `config/`, `include/`, `src/`, `tests/`, and a
Makefile. Shared exposure-list, path, MPI, scheduler, and numerical utilities
live below `include/general/` and `src/general/`; phase code lives below
`process_astrocat`, `process_extcat`, `process_init`, `process_main`,
`process_rearr`, and `process_fd`.

## Top-level phases

The executable invokes enabled phases in this fixed order:

| Phase | CLI switch | Purpose |
|---|---|---|
| `process_astrocat` | `--run-astrocat` | Repartition raw two-column Gaia catalogs into deduplicated one-degree tiles. |
| `process_extcat` | `--run-extcat` | Repartition raw catalogs into the pipeline's standard tile format. |
| `process_init` | `--run-init` | Extract Science images and DQ-mask chips and publish exposure lists. |
| `process_main` | `--run-main` | Run the nine numerical stages. |
| `process_rearr` | `--run-rearr` | Spatially partition `*_all.cat`. |
| `process_fd` | `--run-fd` | Perform the field-distortion shear test. |

`process_astrocat` and `process_extcat` run once, in that order. Other phases
run per dataset; datasets are sequential and the first failure stops the run.

`process_main` stages are selected by prime factors in `PROCESS_stage`:

| Stage | Prime | Work |
|---:|---:|---|
| 1 | 2 | preprocessing and Gaia matching |
| 2 | 3 | astrometry |
| 3 | 5 | source detection and star candidates |
| 4 | 7 | star-candidate power spectra |
| 5 | 11 | PSF modeling |
| 6 | 13 | galaxy power spectra |
| 7 | 17 | Fourier_Quad shear estimators |
| 8 | 19 | exposure diagnostics |
| 9 | 23 | catalog combination and calibration |

The default `223092870` enables all stages. Stage 9 requires Stage 8.

## Build

Use an MPI C++ wrapper with C++17, CFITSIO, FFTW3, Eigen3, LAPACK, and BLAS:

```bash
cd cpp_Lite                    # or cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
```

For alternate installations:

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack \
     EIGEN_INCLUDE=/opt/eigen/include/eigen3 -j4
```

The current Makefiles expose only `all` and `clean`. Validate a build first with
`./Fourier_Quad_Pipe --help`, then run a representative phase or dataset for the
configuration being changed. The repository's pinned container stack uses GCC
12.3.0, OpenMPI 4.1.8, CFITSIO 4.6.4, FFTW3 3.3.11, and Eigen3 3.4.0.

## Configuration and CLI

All compiled defaults and fixed scientific parameters are stored in the selected
variant's `config/` directory. Review them before a run; editing a configuration
header requires rebuilding the executable.

The CLI provides per-run overrides for phase control, I/O directories, catalog
schema fields, and the other values represented by
`ProcessConfig::RuntimeOptions`. It does not override most scientific parameters
in `config/LensingConfig.hpp`. Use `./Fourier_Quad_Pipe --help` for the exhaustive
current option list.

Options accept `--name value` and `--name=value`. Booleans accept
`true/false`, `1/0`, and `on/off`. The first explicit `--dataset`, `--contains`,
or `--extcat-contains` replaces its compiled list; repeats append. One bare
argument is accepted as a legacy `--expo-list` alias.

### Common and data-source-dependent parameters

Review the entries below before a run. “Runtime” means the listed CLI can change
the value without rebuilding; “compile-time” means edit the selected variant and
run `make` again. Do not edit derived dimensions or column indices independently.

| Category | Parameters (current defaults) | How to change | When to change / constraints |
|---|---|---|---|
| Top-level phases | `RUN_PROCESS_ASTROCAT`, `RUN_PROCESS_EXTCAT`, `RUN_PROCESS_INIT`, `RUN_PROCESS_MAIN`, `RUN_PROCESS_REARR`, `RUN_PROCESS_FD` | `config/ProcessConfig.hpp`; runtime `--run-astrocat`, `--run-extcat`, `--run-init`, `--run-main`, `--run-rearr`, `--run-fd` | Select work for this invocation. Standard defaults to `false/false/true/true/true/true`; Lite to `false/false/true/true/false/false`. |
| Science/DQ archives and datasets | `SCIENCE_ROOT`, `DQ_ROOT`, `OUTPUT_ROOT`, `DATASETS`, `CONTAINS` | `config/InitConfig.hpp`; runtime `--science-root`, `--dq-root`, `--output-root`, `--dataset`, `--contains` | Change for another observation archive, basename prefix, discovery token, or output root. Lite requires per-chip DQ masks. |
| Exposure lists and phase outputs | `EXPO_LIST`, `REARR_OUTPUT_DIRECTORY`, `REARR_OUTPUT_BASE_DIRECTORY`, `REARRANGED_EXPO_LIST_FILENAME`, `REARRANGED_EXPO_LIST_DIRECTORY`, `FD_EXPO_LIST`, `FD_OUTPUT_DIRECTORY`, `FD_OUTPUT_BASE_DIRECTORY` | `config/ProcessConfig.hpp`; runtime `--expo-list`, `--rearr-output-dir`, `--rearr-output-base`, `--rearr-list-name`, `--rearr-list-dir`, `--fd-expo-list`, `--fd-output-dir`, `--fd-output-base` | Change for downstream-only runs or different rearrangement/FD list and output locations. |
| Gaia-catalog tiling | `ASTROCAT_INPUT_DIRECTORY`, `ASTROCAT_OUTPUT_DIRECTORY`, `ASTROCAT_ADD_HEADER=true`, `ASTROCAT_EXISTING_POLICY=fail` | `config/AstroCatConfig.hpp`; runtime `--astrocat-input`, `--astrocat-output`, `--astrocat-add-header`, `--astrocat-existing` | Change for a new raw Gaia catalog or rerun policy. The output option controls only `process_astrocat`; it is not checked against or propagated to `ASTROMETRY_CAT`. |
| Gaia catalog layout | `AstroCatType=1` | `config/LensingConfig.hpp`, compile-time | `1` reads legacy large `gaia_*.cat` tiles; `2` accumulates the one-degree `des_y6_*.dat` tiles produced by `process_astrocat`. Set `ASTROMETRY_CAT` separately to the directory consumed by Stage 1 and rebuild after changing the type. |
| External-catalog discovery and parsing | `EXTCAT_INPUT_DIRECTORY`, `EXTCAT_OUTPUT_DIRECTORY` | `config/ExtCatConfig.hpp`; runtime `--extcat-input`, `--extcat-output` | Change the external-catalog file organization. The output directory cannot equal or sit below the input directory. |
| External-catalog schema | `EXTCAT_TOTAL_COLUMNS`, `EXTCAT_INPUT_COLUMNS_ONE_BASED`, `EXTCAT_RA_COLUMN_ONE_BASED`, `EXTCAT_DEC_COLUMN_ONE_BASED`, `EXTCAT_ZP_COLUMN_ONE_BASED` | `config/ExtCatConfig.hpp`; runtime `--extcat-columns`, `--extcat-ra-column`, `--extcat-dec-column`, `--extcat-zp-column` for projection/field locations | Change for another survey or column order. A projection must retain RA, Dec, ZP, and fields consumed by enabled phases; changing total width also requires review of rearrangement and FD indices. |
| Gaia, source-catalog, and calibration paths | `ASTROMETRY_CAT`, `SOURCE_CAT_DEFAULT` (effective `SOURCE_CAT`), `FLAT_PATH`, `PSF_PATH` | `config/LensingConfig.hpp`; `--extcat-output` sets the effective `SOURCE_CAT` at runtime; the others are compile-time | Change for another Gaia tile set, normalized source catalog, flat, or external-PSF source. `--astrocat-output` remains independent of `ASTROMETRY_CAT`. Container paths must match bind destinations. |
| Standard branch selection | `ASTROMETRY_trivial=0`, `include_FLAT=0`, `include_Mask=2`, `ext_cat=1`, `ext_PSF=0`, `PSF_type=1`, `PSF_Ms=0` | `config/LensingConfig.hpp`, compile-time | Only Standard can switch these branches. Lite is fixed to Gaia, no flat, per-chip DQ, external source catalog, frame-star PSF, local polynomial, and no PCA. |
| Image and detector geometry | `npx=3000`, `npy=5000`, `CCD_split=2`, `chipnx=2046`, `chipny=4094`, `pixel_size=0.2628`, `NMAX_CHIP=62`, `NMAX_EXPO=25000` | `config/LensingConfig.hpp`, compile-time | Change for another camera, CCD size, amplifier layout, pixel scale, or batch exposure count; review geometry values as a coupled set. |
| Numerical stages | `PROCESS_stage=223092870` | `config/LensingConfig.hpp`, compile-time | Select the nine main stages by prime factors; Stage 9 (23) requires Stage 8 (19). |
| Detection and pixel thresholds | `saturation_thresh=25000` | `config/LensingConfig.hpp`, compile-time | Recalibrate with representative data after changing the image source. |
| FD catalog layout | `col_flags_*`, `col_cra/cdec`, `col_mag_*`, `col_zp`, `col_expo`, `col_ccd`, and derived `col_*` | `config/FDConfig.hpp`, compile-time with coordinated reader/writer changes | Defaults assume the 18-field DES schema followed by `EXPO_NUM` and `ccD_NUM`. If external width/order changes, review `ExtCatConfig`, the external reader, rearrangement layout, and FD reader together; never change one index alone. |
| FD detector rules | `bad_ccds={2,31,53,61}`, `chip_xmin=50`, `chip_xmax=1990`, `chip_ymin=100`, `chip_ymax=3990` | `config/FDConfig.hpp`, compile-time | Change for another camera, bad-CCD list, or edge-mask policy; keep `n_bad_ccds` synchronized. |

See [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) for every individual
parameter's Standard/Lite default, legal values, CLI mapping, and rebuild requirement.
For coupled changes, preserve a baseline configuration and validate first on the
smallest representative dataset.

## Run examples

Main-only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-main true \
  --expo-list /data/work/expo_gband.list
```

Initializer plus main:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init true --run-main true \
  --science-root /data/archive/science \
  --dq-root /data/archive/dqmask \
  --output-root /data/work \
  --dataset g2019:c4d_19 --existing resume
```

Repartition raw Gaia catalogs only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-astrocat true --run-extcat false --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --astrocat-input /data/raw_gaia --astrocat-output /data/gaia/tiles \
  --astrocat-add-header true --astrocat-existing fail
```

The `--astrocat-output` path above controls only this phase. Set
`LensingConfig::ASTROMETRY_CAT` separately when `process_main` should consume
that directory.

Normalize an External source catalog for this pipeline only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

At least one phase must be enabled. Successful initialization supplies the
generated absolute `expo_<target>.list` to later phases. Multiple datasets run
sequentially and independently.

## Inputs and outputs

Prepare the four user input classes—Science images, Gaia catalog, External
source catalog, and configuration-dependent DQ masks—according to the
[top-level input data requirements](README.md#input-data-requirements). That
section is the single minimum-schema contract; this guide describes the C++
runtime layout and products.

An exposure list contains one chip-list path per nonblank record; a trailing
legacy chip count is accepted. Initialization reads archives in place and
creates `science/`, `dqmask/`, `stamps/`, and `result/` below each dataset,
plus top-level exposure/fits lists and a manifest.

Principal products are:

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

`_all.cat` is a per-exposure shear catalog. By default, each row contains the
external source catalog fields, the original 1-based `EXPO_NUM`, one CCD number,
and 25 pipeline fields (45 columns total). Regenerate Stage 9, rearrangement,
and FD products after this schema change; legacy 44-column rows are incompatible.
`subcat_*.cat` spatially repartitions those rows by Dec and RA; repeated
measurements at identical coordinates are adjacent, which simplifies
deduplication. `FD_test_comb.dat` is the field-distortion test table written by
`process_fd` and is used to calibrate shear measurements.

Common errors include:

- Do not enable Stage 9 without Stage 8.
- The external-catalog output cannot equal or be below its input directory.
- Explicit projection must preserve the fields consumed by every enabled stage.
- Lite cannot enable branches that were removed from its source.
- Use container paths in CLI arguments executed inside Docker/Apptainer.

For settings and the output schema, see
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md). For containers, see
[cpp_docker/README.md](cpp_docker/README.md) and the
[Slurm runner](cpp_docker/runner/README.md).
