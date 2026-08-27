# Fourier_Quad C++ Pipeline Guide

This guide describes the C++ programs in this repository. They use CLI
overrides and compiled defaults; unlike the separate `Fourier_Quad_Cpp`
repository, this version has no INI configuration layer.

> 中文版：[CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)

## Variants and layout

[`cpp_Standard`](cpp_Standard/) retains optional flat, mask, identity
astrometry, external/hybrid PSF, and PCA branches. [`cpp_Lite`](cpp_Lite/)
physically removes those alternatives and keeps Gaia astrometry, per-chip DQ
masks, external sources, deblending, local-polynomial PSF, and no PCA.

Each variant contains `main.cpp`, `config/`, `include/`, `src/`, `tests/`, and a
Makefile. Shared exposure-list, path, MPI, scheduler, and numerical utilities
live below `include/general/` and `src/general/`; phase code lives below
`process_extcat`, `process_init`, `process_main`, `process_rearr`, and
`process_fd`.

## Top-level phases

The executable invokes enabled phases in this fixed order:

| Phase | CLI switch | Purpose |
|---|---|---|
| `process_extcat` | `--run-extcat` | Repartition raw text catalogs into sky tiles. |
| `process_init` | `--run-init` | Extract Science/DQ chips and publish exposure lists. |
| `process_main` | `--run-main` | Run the nine numerical stages. |
| `process_rearr` | `--run-rearr` | Spatially partition `*_all.cat`. |
| `process_fd` | `--run-fd` | Perform the field-distortion shear test. |

`process_extcat` runs once. Other phases run per dataset; datasets are
sequential and the first failure stops the run.

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

Current Make targets are `all`, `clean`, and
`test-general-infrastructure`.

## Configuration and CLI

CLI values override fields represented by `ProcessConfig::RuntimeOptions`.
Most scientific values in `config/LensingConfig.hpp`, rearrangement constants,
and FD statistics remain compile-time and require rebuilding.

Options accept `--name value` and `--name=value`. Booleans accept
`true/false`, `1/0`, and `on/off`. The first explicit `--dataset`, `--contains`,
or `--extcat-contains` replaces its compiled list; repeats append. One bare
argument is accepted as a legacy `--expo-list` alias.

Use `./Fourier_Quad_Pipe --help` for the exhaustive current list. The main
groups are:

- five `--run-*` phase switches;
- `--extcat-*` input, schema, chunking, and publication options;
- `--science-root`, `--dq-root`, `--output-root`, repeatable datasets/archive
  filters, and the initializer existing-output policy;
- `--expo-list`, `--rearr-*`, and `--fd-*` downstream paths.

`--extcat-output` also changes the effective external source-catalog path for
that invocation. Other lensing paths and branch choices remain compiled.

## Run examples

Main-only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

Initializer plus main:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init true --run-main true --run-rearr false --run-fd false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dqmask \
  --output-root /data/work \
  --dataset g2019:c4d_19 --existing resume
```

External-catalog-only:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

At least one phase must be enabled. Successful initialization supplies the
generated absolute `expo_<target>.list` to later phases. Multiple datasets run
sequentially.

## Inputs and outputs

An exposure list contains one chip-list path per nonblank record; a trailing
legacy chip count is accepted. Initialization reads archives in place and
creates `science/`, `dqmask/`, `stamps/`, and `result/` below each dataset,
plus top-level exposure/fits lists and a manifest.

Principal products are:

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<rearr-output-dir>/catalog_summary.txt
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

Stage 7 writes 24 fields through WCS parity. Stage 9 appends exposure
chi-square, so the default final row has 18 external fields + one CCD number +
25 pipeline fields = 44 fields. Explicit external projection changes the
external prefix width. RA, Dec, and photo-z must remain available to the
pipeline.

## Common errors

- Do not enable Stage 9 without Stage 8.
- The external-catalog output cannot equal or be below its input directory.
- Explicit projection must contain raw RA, Dec, and photo-z fields.
- Lite cannot enable branches that were removed from its source.
- Use container paths in CLI arguments executed inside Docker/Apptainer.

For settings and the output schema, see
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md). For containers, see
[cpp_docker/README.md](cpp_docker/README.md) and the
[Slurm runner](cpp_docker/runner/README.md).
