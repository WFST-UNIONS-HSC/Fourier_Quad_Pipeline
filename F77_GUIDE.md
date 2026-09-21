# Fourier_Quad Fortran Pipeline Guide

The `f77` and `f77_Lite` directories contain the legacy MPI pipeline. This
program is configured entirely at build time and accepts one positional
exposure-list path.

> 中文版：[F77_GUIDE_CN.md](F77_GUIDE_CN.md)

## Choose a variant

- `f77` (F77 Standard): full branch set, including optional PCA/multi-scale PSF storage in
  `00_psf_module.f`.
- `f77_Lite` (F77 Lite): fixed production path with alternate astrometry, flat, mask,
  source, PSF, deblending, hybrid, and PCA branches removed.

Both variants build `Fourier_Quad_Pipe`.

## Numerical stages

`PROCESS_stage` in `para.inc` is a product of stage primes:

| Stage | Prime | Fortran entry | Work |
|---:|---:|---|---|
| 1 | 2 | `pre_process` | background/noise preprocessing |
| 2 | 3 | `proc_astrometry` | Gaia astrometry |
| 3 | 5 | `proc_source` | source detection and star candidates |
| 4 | 7 | `proc_FourierT_st1` | star-candidate power spectra |
| 5 | 11 | `proc_PSF` | PSF modeling |
| 6 | 13 | `proc_FourierT_st2` | galaxy power spectra |
| 7 | 17 | `proc_shear` | Fourier_Quad shear estimators |
| 8 | 19 | `proc_info` | exposure diagnostics |
| 9 | 23 | `proc_comb` | catalog combination and calibration |

The default product `223092870` enables all stages. Keep Stage 8 enabled when
running Stage 9 so the combined catalog has valid exposure diagnostics.

## Configure

Prepare Science images, Gaia catalog, External source catalog, and
configuration-dependent DQ masks according to the
[top-level input data requirements](README.md#input-data-requirements). F77 Lite
uses its frozen mask path, so do not omit DQ masks unless the selected source and
configuration have been verified not to read them.

All settings require rebuilding:

| File | Purpose |
|---|---|
| `para.inc` | stage selector, catalog/calibration paths, stamp geometry, branch controls, thresholds, limits, and catalog indices |
| `cust_para.inc` | CCD geometry and Standard PCA settings |
| `sig_para.inc` | robust mode-bar noise-plane estimator |

At minimum, review `PROCESS_stage`, `ASTROMETRY_CAT`, `SOURCE_CAT`, and any
active `FLAT_PATH`. In a container, these strings must be container paths that
match the bind destinations.

`SOURCE_CAT_TILE_PREFIX` controls the external-catalog tile basename prefix
(the default is `extern_`). `strl` is 512 in both variants; longer exposure-list,
chip-list, input, or generated product paths are rejected instead of truncated.

Lite documents its frozen branch behavior at the top of `para.inc`; the
deleted alternatives cannot be re-enabled by adding a parameter.

## Build

Required tools are `mpif77`, CFITSIO, LAPACK, and BLAS; Fourier routines are
compiled from the included `FFTPACK.f`, so the current Makefile does not require
an external FFTW library. See the top-level F77 prerequisite table for regular
Linux and HPC requirements.

```bash
cd f77                         # or f77_Lite
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

If CFITSIO is not named `libcfitsio.so` in that directory, pass its full path
as `CFITSIO_LIB`. The Makefiles expose `all` and `clean`. Their compiled default
library directories are site-specific, so portable and container builds
should pass the overrides explicitly.

## Exposure list and run

Each exposure-list record must contain the chip-list path and chip count:

```text
/data/work/stamps/123456.list 60
/data/work/stamps/123457.list 59
```

Blank records are ignored in both the exposure list and each chip list. Science
paths in a chip list must follow `<dataset>/science/<exposure>/<chip>.fits` so
the pipeline can derive `<dataset>` by removing three path components. Exposure
and chip counts are checked against `NMAX_EXPO` and `NMAX_CHIP`.

Run with one positional argument:

```bash
mpirun -np 8 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

There is no `--help`, `--run-*`, or `--config` interface in the Fortran
program. Select stages in `para.inc` and rebuild.

## Outputs

Stages write their intermediate products below the dataset tree referenced by
the chip lists. Standard and Lite now use the same case-sensitive product
layout as `cpp_Standard` and `cpp_Lite`:

| Product family | Location |
|---|---|
| DQ mask | `dqmask/<exposure>/<exposure>_<ccd>.fits` |
| Normalized images and astrometry | `stamps/Norm/<exposure>/`, `astrometry/dat_Astro/<exposure>/`, `astrometry/Head/`, `astrometry/dat_Chk/` |
| Source and star intermediates | product-specific `stamps/cat_Orig`, `stamps/dat_*`, and `stamps/fits_*` directories, with an `<exposure>/` subdirectory for chip products |
| Exposure PSF products | `stamps/dat_StarInfo`, `stamps/fits_StarP`, `stamps/fits_PsfSrc`, `stamps/dat_StarComp`, and `stamps/dat_Rescale` |
| Standard PCA products | `stamps/dat_StarXY/<exposure>/`, `stamps/fits_PsfResi/<exposure>/`, `stamps/dat_Pcs/`, and `stamps/dat_StarCompV2/` |
| Final catalog | `result/<exposure>_all.cat` |

The program does not create directories. Initialize the dataset tree with the
C++ initializer, or create the complete product hierarchy before running F77.
Stage 8 also writes the parent-list summary `expo_info.dat` beside the exposure
list. Every Stage 9 row contains one-based `EXPO_NUM`, physical `ccD_NUM`, the
24-field shear record, and `Chi2`; external-catalog columns precede these fields
when enabled.

Do not modify original archives or catalogs in place. Use a writable processing
tree and ensure every MPI rank can see the same absolute paths.

## Containers and Slurm

For a reproducible local toolchain, see
[f77_docker/README.md](f77_docker/README.md). For non-root Slurm deployment,
see [f77_docker/runner/README.md](f77_docker/runner/README.md). The image does
not include the source or data; both are bind-mounted.
