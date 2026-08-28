# Fourier_Quad Pipeline

> 中文版：[README_CN.md](README_CN.md)

## Overview

Modern C++17 and Legacy F77 implementation of the Fourier_Quad weak-lensing shear measurement pipeline. The
repository provides current C++17 and legacy Fortran implementations, each in Standard and Lite variants.

## Pipeline variants

| Pipeline | Source directory | Purpose |
|---|---|---|
| C++ Lite | [`cpp_Lite`](cpp_Lite/) | C++17 production path with unused alternate branches removed. |
| C++ Standard | [`cpp_Standard`](cpp_Standard/) | Full C++17 branch set, including optional scientific paths. |
| F77 Lite | [`f77_Lite`](f77_Lite/) | Reduced legacy Fortran production path. |
| F77 Standard | [`f77`](f77/) | Full legacy Fortran branch set. |

All four variants build an executable named `Fourier_Quad_Pipe`.

## Quick Start

### 1. Download a Release source package

Open [GitHub Releases](https://github.com/WFST-UNIONS-HSC/Fourier_Quad_Pipeline/releases)
and download only the source package for the Pipeline you plan to run:

| Pipeline | Release source package |
|---|---|
| C++ Lite | `cpp_Lite.zip` |
| C++ Standard | `cpp_Standard.zip` |
| F77 Lite | `f77_Lite.zip` |
| F77 Standard | `f77.zip` |

Use a fixed Release so the software version used for an analysis can be
recorded and reproduced.

### 2. Prepare the environment

Install the dependencies in the relevant [C++](#c-pipeline) or
[F77](#f77-pipeline) environment table below.

### 3. Prepare input data

Prepare the four input classes described in
[Input data requirements](#input-data-requirements). The first three are
required; DQ masks depend on the selected Pipeline configuration.

### 4. Configure the Pipeline

- C++ users: review the selected source package's `config/*.hpp`. The complete
  Standard/Lite comparison is in
  [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md).
- F77 users: edit `para.inc`, `cust_para.inc`, and `sig_para.inc` as described
  in the [F77 guide](F77_GUIDE.md).

### 5. Build and run

From the extracted source-package directory, the shortest entry points are:

```bash
# C++ Lite or C++ Standard
make -j4
./Fourier_Quad_Pipe --help

# F77 Lite or F77 Standard
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
mpirun -np 4 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

See the [C++ guide](CPP_GUIDE.md) or [F77 guide](F77_GUIDE.md) for complete
arguments, phase selection, and run examples.

## Environment prerequisites

### C++ Pipeline

| Category | Prerequisite | Regular Linux | HPC / Slurm | Notes |
|---|---|---:|---:|---|
| Operating system | 64-bit Linux | Required | Required | A modern Linux distribution is recommended. |
| C++ compiler | MPI C++ wrapper with C++17 support | Required | Required | The current container toolchain uses GCC 12.3.0. |
| MPI | OpenMPI or a compatible MPI implementation | Required | Required | The current container uses OpenMPI 4.1.8. |
| CFITSIO | CFITSIO development library | Required | Required | FITS/FZ I/O; container version 4.6.4. |
| FFTW3 | Double- and single-precision FFTW3 libraries | Required | Required | Fourier transforms; container version 3.3.11. |
| Eigen3 | Eigen3 headers | Required | Required | C++ numerical operations; container version 3.4.0. |
| BLAS / LAPACK | BLAS and LAPACK libraries | Required | Required | Linear algebra. |
| Build tools | `make`, shell, standard GNU tools | Required | Required | Builds and helper scripts. |
| Shared filesystem | All ranks/nodes see the same inputs and outputs | No | Required | Required for multi-node jobs. |
| Slurm | Site scheduler and `srun` | No | Required for Slurm jobs | Batch allocation and launch. |
| PMI2-compatible launch | Slurm `pmi2` support and compatible MPI | No | Required for the repository HPC container path | Used for direct Slurm MPI launch. |
| Apptainer / Singularity | Rootless container runtime | No | Required when using the repository HPC container path | Must be available on compute nodes. |
| Writable processing directory | Shared output/work directory writable by all ranks | Required | Required | Source data and outputs should not be mixed. |

### F77 Pipeline

| Category | Prerequisite | Regular Linux | HPC / Slurm | Notes |
|---|---|---:|---:|---|
| Operating system | 64-bit Linux | Required | Required | Target platform for the F77 Pipeline. |
| Fortran compiler | Compatible GNU Fortran toolchain with `mpif77` | Required | Required | The pinned container uses GNU Fortran 4.8.5. |
| MPI | MPICH or a compatible Fortran MPI implementation | Required | Required | The pinned container uses MPICH 4.1.2. |
| CFITSIO | CFITSIO library | Required | Required | FITS/FZ I/O; container version 4.3.1. |
| Fourier routines | Repository-provided `FFTPACK.f` | Included | Included | The current F77 Makefile does not require an external FFTW library. |
| BLAS / LAPACK | BLAS and LAPACK libraries | Required | Required | The container supplies LAPACK 3.8.0 and its BLAS. |
| Build tools | `make`, shell, standard GNU tools | Required | Required | Builds and run scripts. |
| Shared filesystem | All tasks/nodes see the same inputs and outputs | No | Required | Required for multi-node jobs. |
| Slurm | Site scheduler | No | Required for Slurm jobs | Batch allocation and launch. |
| MPI-compatible launch | A launcher compatible with the selected container/host MPI mode | No | Required | The repository runner supports its documented MPICH and PMI-compatible modes. |
| Apptainer / Singularity | Rootless container runtime | No | Required when using the repository HPC container path | Must be available on compute nodes. |
| Writable processing directory | Shared output/work directory writable by all tasks | Required | Required | Do not write products into original archives. |

## Input data requirements

| Input | Required | Purpose | Minimum requirement |
|---|---:|---|---|
| Science images | Yes | Exposures on which the Pipeline performs source detection, shape measurement, and weak-lensing processing. | Supported Science FITS/FZ data whose file organization matches the configured exposure and CCD recognition rules. |
| Gaia catalog | Yes | High-precision sky-coordinate reference for source matching and astrometric calibration. | Covers the Science-image footprint; each file has one header line, followed by rows whose two fields are numeric `ra` and `dec`. |
| External source catalog | Yes | Supplies sky positions, `zp`, and photometry for external-source matching and downstream processing. | Contains `ra`, `dec`, `zp`, and a magnitude in at least one observed band. |
| DQ masks | Configuration-dependent (optional input class) | Marks bad, saturated, defective, or otherwise invalid pixels. | May be omitted only when the selected configuration does not read DQ masks. |

### Science images

Science images are the primary scientific exposures, not calibration catalogs or
an output directory. They must use a FITS/FZ format supported by the selected
variant, be readable through its FITS logic, and follow the exposure/CCD naming
and list conventions described in the detailed guide.

### Gaia catalog

The Gaia catalog supplies accurate RA/Dec reference positions for object
matching and astrometric calibration. It must cover the actual Science-image
footprint and be stored directly under the configured `ASTROMETRY_CAT` directory.
The reader skips the first line as a header, then reads the first two numeric
fields of each remaining row as RA and Dec. Rows may be comma- or
whitespace-separated; additional fields are ignored.

**Filename convention:**

- `|Dec| < 80°`: `gaia_<p|m><D>_<RR>.cat`, where
  `D = floor(|Dec| / 10) + 1` (1-8) and `RR = floor(RA / 10)` (00-35,
  zero-padded).
- `|Dec| >= 80°`: `gaia_<p|m>9.cat`, without an RA suffix.
- `p` denotes nonnegative Dec; `m` denotes negative Dec.
- Each file must contain one header line.

> Examples:
> 1. gaia_p1_00.cat covers `0° <= RA < 10°` and `0° <= Dec < 10°`
> 2. gaia_m3_12.cat covers `120° <= RA < 130°` and `-30° < Dec <= -20°` 
> 3. gaia_p9.cat covers `0° <= RA < 360°` and `80° <= Dec <= 90°`
*To ensure that stars can still be selected for exposures located at the edges of the 10°×10° grid, it is recommended to expand the upper and lower Dec coverage limits of a single star catalog by 2° based on the aforementioned limits, and expand the upper and lower RA limits by 2°, 4°, and 6° within the 0°, 30°, and 60° ranges, respectively.*
### External source catalog

The minimum schema is intentionally survey- and band-independent:

| Field | Meaning |
|---|---|
| `ra` | Right Ascension. |
| `dec` | Declination. |
| `zp` | The catalog `zp` quantity consumed by the selected Pipeline configuration. |
| One observed-band magnitude | A magnitude in any one band used by the selected analysis. |

Additional colors, redshifts, object classes, shapes, and flags may be retained,
but they are not part of this minimum input contract. Configure actual column
positions, names, delimiter, header handling, and projection in
`config/ExtCatConfig.hpp` or with its documented CLI overrides.

**Filename convention:**

- 1° × 1° tile: `des_y6_RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`.
- RA boundaries use three digits. Dec boundaries use `p` or `m` plus a two-digit
  absolute value. Each upper boundary is one degree above its lower boundary.
- Each file must contain one header line.

> Example:
> `des_y6_RA_123_124_Dec_m05_m04.dat` covers `123° <= RA < 124°` and
> `-5° <= Dec < -4°`.

### DQ masks (optional)

DQ masks identify pixels that must not participate in scientific measurements,
including bad pixels, saturation, detector defects, and other invalid regions.
They are optional only when the chosen configuration disables DQ access. Standard
users can select the relevant mask mode with `include_Mask`; current C++ Lite is
fixed to per-chip DQ masks, so Lite runs must provide them. If DQ masks are
omitted, verify that no configured path or active branch still reads them.

## C++ configuration reference

Use [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) to compare every
configuration header, Standard/Lite default, CLI override, and rebuild rule.

## Detailed guides

- [C++ Pipeline guide](CPP_GUIDE.md)
- [F77 Pipeline guide](F77_GUIDE.md)
- [External source catalog tooling](gen_src_cat/README.md)

## License

Repository-authored code is distributed under the [MIT License](LICENSE).
Container dependencies retain their upstream licenses; consult each container
directory's third-party notices.
