# Fourier_Quad_Pipeline

A weak-lensing shear measurement pipeline based on the **Fourier\_Quad** method,
provided in both **Fortran 77** and **C++** implementations, together with
ready-to-use **Docker** build environments and **HPC** (Slurm/Apptainer) runner
scripts.

> **中文文档**：请参阅 [README_CN.md](README_CN.md)

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Pipeline Stages](#pipeline-stages)
- [Configuration](#configuration)
- [Building from Source](#building-from-source)
- [Docker Environments](#docker-environments)
- [HPC Deployment](#hpc-deployment)
- [CI/CD](#cicd)
- [Attribution](#attribution)
- [License](#license)

---

## Overview

The Fourier\_Quad method measures gravitational shear by computing fourth-order
moments (the "quad" term) in Fourier space. This repository contains a complete,
MPI-parallel pipeline that processes astronomical CCD images end-to-end: from
pre-processing and astrometric calibration through source detection, PSF
modeling, and shear estimation, to final catalog combination.

Two independent codebases are maintained:

| Codebase | Language | Compiler | MPI | Key libraries |
|---|---|---|---|---|
| `f77` / `f77_Lite` | Fortran 77 (+ Fortran 90 module) | GCC 4.8.5 (gfortran) | MPICH 4.1.2 | CFITSIO 4.3.1, LAPACK 3.8.0, FFTPACK |
| `cpp_Standard` / `cpp_Lite` | C++17 | G++ 12.3.0 | OpenMPI 4.1.8 | CFITSIO 4.6.4, FFTW 3.3.11, Eigen 3.4.0, LAPACK 3.11.0 |

Each codebase has two variants:

- **Standard** (`f77`, `cpp_Standard`): full feature set, including multi-scale /
  PCA PSF reconstruction (`PSFRecons` / `proc_psfreconsV2.f`).
- **Lite** (`f77_Lite`, `cpp_Lite`): frozen-branch simplified version. Eight
  compile-time switches are fixed to their production values and all dead-code
  branches are physically removed. See `cpp_Lite/REFACTOR_NOTES.md` for details.

---

## Repository Structure

```
Fourier_Quad_Pipeline/
├── f77/                  Full Fortran 77 pipeline source
├── f77_Lite/             Simplified Fortran 77 pipeline (frozen branches)
├── cpp_Standard/         Full C++17 pipeline source
├── cpp_Lite/             Simplified C++17 pipeline (frozen branches)
├── f77_docker/           Docker build environment for f77 pipeline
├── cpp_docker/           Docker build environment for C++ pipeline
├── .github/workflows/    CI/CD workflows (GHCR images + GitHub releases)
├── README.md             This file (English)
├── README_CN.md          Chinese documentation
├── LICENSE               MIT License
└── .gitignore
```

### Source directories

#### `f77/` — Full Fortran 77 pipeline

| File | Description |
|---|---|
| `main.f` | Main program entry point. Initializes MPI, reads the exposure list, and dispatches pipeline stages. |
| `para.inc` | Master parameter file. Controls image dimensions, stage selection (`PROCESS_stage`), catalog paths, PSF order, stamp size, thresholds, and catalog column indices. |
| `cust_para.inc` | Custom parameters for CCD geometry and PCA PSF decomposition. |
| `sig_para.inc` | F6 mode-bar noise-plane estimator parameters. |
| `pre_process.f` | **Stage 1**: flat-field, mask handling, background/noise estimation. |
| `proc_astrometry.f` | **Stage 2**: astrometric calibration using Gaia reference catalog. |
| `proc_source.f` | **Stage 3**: source detection and stamp extraction. |
| `proc_FFT_st1.f` | **Stage 4**: first-stage Fourier transform of galaxy stamps. |
| `proc_PSF.f` | **Stage 5**: PSF modeling from stellar stamps (local polynomial fit). |
| `proc_psfreconsV2.f` | PSF PCA reconstruction (multi-scale, `PSF_Ms=1` only). |
| `00_psf_module.f` | Fortran 90 module for global PSF PCA storage. |
| `proc_FFT_st2.f` | **Stage 6**: second-stage Fourier transform. |
| `proc_shear.f` | **Stage 7**: Fourier\_Quad shear estimation. |
| `proc_info.f` | **Stage 8**: per-exposure statistics collection. |
| `proc_combine_shear_catalog.f` | **Stage 9**: shear catalog combination and calibration. |
| `astrometry_calib.f` | Astrometric calibration utility routines. |
| `ex_star.f` | Star extraction and classification. |
| `rw_fits_image.f` | FITS image read/write routines (via CFITSIO). |
| `universal.f` | General-purpose utility functions. |
| `univ_imag_proc.f` | Image processing utilities. |
| `press.f` | Numerical Recipes routines (sorting, statistics, interpolation). |
| `mpi_routines.f` | MPI distribution helper (`mpi_distribute`). |
| `FFTPACK.f` | FFTPACK 5.1 FFT library. |
| `Makefile` | Build file. Uses `mpif77`, links against LAPACK, BLAS, and CFITSIO. |

#### `f77_Lite/` — Simplified Fortran 77 pipeline

Identical file set to `f77/` except `00_psf_module.f` is absent (PCA PSF
reconstruction is removed). All eight compile-time switches are frozen to
production values and dead-code branches are physically removed.

#### `cpp_Standard/` — Full C++17 pipeline

| File | Description |
|---|---|
| `main.cpp` | Main program entry point. Mirrors `f77/main.f` in C++. |
| `LensingConfig.hpp` | Configuration constants (equivalent to `para.inc` + `cust_para.inc` + `sig_para.inc`). |
| `PreProcess.cpp/.hpp` | **Stage 1**: pre-processing. |
| `Astrometry.cpp/.hpp` | **Stage 2**: astrometric calibration. |
| `SourceExtractor.cpp/.hpp` | **Stage 3**: source detection and extraction. |
| `FourierTransformSt1.cpp/.hpp` | **Stage 4**: first-stage Fourier transform. |
| `PSFModel.cpp/.hpp` | **Stage 5**: PSF modeling. |
| `PSFRecons.cpp/.hpp` | PSF PCA reconstruction (`PSF_Ms=1` only). |
| `FourierTransformSt2.cpp/.hpp` | **Stage 6**: second-stage Fourier transform. |
| `ShearMeasurement.cpp/.hpp` | **Stage 7**: Fourier\_Quad shear estimation. |
| `ExposureInfo.cpp/.hpp` | **Stage 8**: per-exposure statistics. |
| `CatalogCombiner.cpp/.hpp` | **Stage 9**: catalog combination and calibration. |
| `FitsIO.cpp/.hpp` | FITS I/O routines. |
| `LinearSolve.cpp/.hpp` | Linear algebra utilities. |
| `UniversalUtils.cpp/.hpp` | General-purpose utilities. |
| `ImageProcessing.cpp/.hpp` | Image processing utilities. |
| `NumericalRecipes.cpp/.hpp` | Numerical Recipes port (RNG, sorting, interpolation). |
| `MPIScheduler.cpp/.hpp` | MPI initialization and task distribution. |
| `ExStar.cpp/.hpp` | Star extraction and classification. |
| `Makefile` | Build file. Uses `mpicxx`, C++17, links against CFITSIO, FFTW, LAPACK. |

#### `cpp_Lite/` — Simplified C++17 pipeline

Same file set as `cpp_Standard/` except `PSFRecons.cpp/.hpp` is absent.
See `cpp_Lite/REFACTOR_NOTES.md` for the detailed change log.

### Docker directories

#### `f77_docker/`

Builds a development image matching the pilogin cluster toolchain:

| Component | Version |
|---|---|
| Base image | Rocky Linux 8.10 (glibc 2.28) |
| GCC / G++ / GFortran | 4.8.5 (compiled from source) |
| MPICH | 4.1.2 (ch4:ofi device, Slurm launcher) |
| CFITSIO | 4.3.1 |
| LAPACK + reference BLAS | 3.8.0 |

Key files: `Dockerfile`, `compose.yaml`, `.env.example`, `scripts/verify-image.sh`,
`scripts/check-public-repo.sh`, `runner/` (HPC deployment scripts), `config/`,
`patches/`, `checksums.sha256`, `SOURCES.md`, `THIRD_PARTY_NOTICES.md`.

#### `cpp_docker/`

Builds a portable HPC runtime image:

| Component | Version |
|---|---|
| Base image | Rocky Linux 8.10 |
| G++ | 12.3.0 (conda-forge) |
| OpenMPI | 4.1.8 (Slurm PMI2 direct-launch) |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / OpenBLAS | 3.11.0 / 0.3.33 |

Key files: `Dockerfile`, `compose.yaml`, `pixi.toml`, `.env.example`,
`scripts/verify-image.sh`, `scripts/check-public-repo.sh`, `runner/` (HPC
deployment scripts), `SOURCES.md`, `THIRD_PARTY_NOTICES.md`.

---

## Pipeline Stages

The pipeline consists of 9 stages. Stage execution is controlled by the
`PROCESS_stage` parameter (defined in `para.inc` / `LensingConfig.hpp`), which is
the product of prime factors. A stage runs when `PROCESS_stage` is divisible by
its prime. The default value `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` enables
all stages.

| Stage | Prime | Function | Description |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | Read FITS images, apply flat-field and mask corrections, estimate background noise (F6 mode-bar estimator). |
| 2 | 3 | `proc_astrometry` / `Astrometry` | Astrometric calibration using Gaia reference catalog; WCS fitting. |
| 3 | 5 | `proc_source` / `SourceExtractor` | Source detection, deblending, and stamp extraction. |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | First-stage Fourier transform of galaxy stamps. |
| 5 | 11 | `proc_PSF` / `PSFModel` | PSF modeling from stellar stamps via local polynomial fitting. Optional PCA reconstruction (`PSF_Ms=1`). |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | Second-stage Fourier transform. |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | Fourier\_Quad shear estimation from fourth-order Fourier moments. |
| 8 | 19 | `proc_info` / `ExposureInfo` | Collect per-exposure statistics (PSF FWHM, star count, etc.). |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | Combine shear catalogs across exposures and apply calibration corrections. |

To disable a stage, divide `PROCESS_stage` by its prime factor. For example,
setting `PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19` (omit 23) skips
catalog combination.

---

## Configuration

### Fortran 77 (`para.inc`, `cust_para.inc`, `sig_para.inc`)

All compile-time parameters live in three include files:

- **`para.inc`** — Image dimensions (`npx`, `npy`), stage control
  (`PROCESS_stage`), catalog paths, PSF parameters, stamp size, detection
  thresholds, and catalog column indices.
- **`cust_para.inc`** — CCD geometry (`chipnx`, `chipny`, `Camera_ccd_num`) and
  PCA PSF decomposition parameters.
- **`sig_para.inc`** — F6 mode-bar noise-plane estimator parameters.

> Paths such as `ASTROMETRY_CAT`, `SOURCE_CAT`, and `FLAT_PATH` must match the
> container mount paths defined in the Docker `.env` files.

### C++ (`LensingConfig.hpp`)

All parameters from the three Fortran include files are consolidated into a
single `LensingConfig` namespace in `LensingConfig.hpp`. Catalog column indices
are shifted to 0-based for C++.

---

## Building from Source

### Fortran 77

Prerequisites: `gfortran`, `mpif77` (MPICH), CFITSIO, LAPACK/BLAS.

```bash
cd f77          # or f77_Lite
# Edit para.inc to set your catalog paths and parameters
make
# Executable: ./Fourier_Quad_Pipe
```

The Makefile accepts override variables:

```bash
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

### C++

Prerequisites: `g++` (C++17), `mpicxx` (OpenMPI or MPICH), CFITSIO, FFTW3,
LAPACK/BLAS, Eigen3.

```bash
cd cpp_Standard   # or cpp_Lite
# Optionally edit LensingConfig.hpp for your catalog paths
make -j4
# Executable: ./Fourier_Quad_Main
```

The `cpp_Standard` Makefile supports an optional `STACK_PREFIX`:

```bash
make STACK_PREFIX=/opt/cppstack -j4
```

### Running the pipeline

```bash
mpirun -np <N> ./Fourier_Quad_Pipe <EXPO_LIST>    # Fortran
mpirun -np <N> ./Fourier_Quad_Main <EXPO_LIST>    # C++
```

`EXPO_LIST` is a text file listing exposure names and chip counts. Each MPI rank
is assigned a subset of exposures to process.

---

## Docker Environments

The Docker environments provide a reproducible build toolchain without requiring
manual installation of compilers and scientific libraries.

### Quick start (f77)

```bash
cd f77_docker
cp .env.example .env
# Edit .env: set F77_SOURCE_HOST, catalog paths, and PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-F77
```

Inside the container:

```bash
cd /workspace/f77
make
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```

### Quick start (C++)

```bash
cd cpp_docker
cp .env.example .env
# Edit .env: set CPP_SOURCE_HOST, catalog paths, and PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-CPP
```

Inside the container:

```bash
cd /workspace/cpp_Standard
make -j4
mpirun -np 4 ./Fourier_Quad_Main /data/DataProcess/expo_list.list
```

### Verifying an image

Each Docker directory includes a verification script:

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
bash cpp_docker/scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

For detailed Docker environment documentation, see:
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)
- [`cpp_docker/README.md`](cpp_docker/README.md) / [`cpp_docker/README-CN.md`](cpp_docker/README-CN.md)

---

## HPC Deployment

Both Docker environments include `runner/` directories with Slurm/Apptainer
deployment scripts. The typical workflow is:

1. Pull the GHCR image and convert to a SIF:
   ```bash
   bash f77_docker/runner/pull-sif.sh    # or cpp_docker/runner/pull-sif.sh
   ```

2. Configure the environment:
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # Edit paths for your cluster
   ```

3. Audit MPI compatibility (read-only):
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. Submit the pipeline:
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

For detailed HPC runner documentation, see:
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)
- [`cpp_docker/runner/README.md`](cpp_docker/runner/README.md) / [`cpp_docker/runner/README-CN.md`](cpp_docker/runner/README-CN.md)

---

## CI/CD

Two GitHub Actions workflows automate image publication and releases.

### Docker images → GHCR

Builds and publishes two images to the GitHub Container Registry:

| Image | GHCR path | Tags |
|---|---|---|
| f77 | `ghcr.io/<owner>/<repo>/f77pipeline` | `gnu4.8.5`, `<version>`, `latest`, `<sha>` |
| cpp | `ghcr.io/<owner>/<repo>/cpppipeline` | `gxx12.3-openmpi4.1.8-pmi2`, `<version>`, `latest`, `<sha>` |

Pull a published image:

```bash
docker pull ghcr.io/syoong-s/fourier_quad_pipeline/f77pipeline:latest
docker pull ghcr.io/syoong-s/fourier_quad_pipeline/cpppipeline:latest
```

## Attribution
This project implements the method proposed by Prof. Zhang in:

- Zhang, J. (2007). Measuring the cosmic shear in Fourier space: Measuring the cosmic shear in Fourier space. Monthly Notices of the Royal Astronomical Society, 383(1), 113–118. https://doi.org/10.1111/j.1365-2966.2007.12585.x
- Zhang, J., Luo, W., & Foucaud, S. (2015). Accurate shear measurement with faint sources. Journal of Cosmology and Astroparticle Physics, 2015(01), 024–024. https://doi.org/10.1088/1475-7516/2015/01/024
- Zhang, J., Zhang, P., & Luo, W. (2017). APPROACHING THE CRAMÉR–RAO BOUND IN WEAK LENSING WITH PDF SYMMETRIZATION. The Astrophysical Journal, 834(1), 8. https://doi.org/10.3847/1538-4357/834/1/8


## License

Repository-authored files are licensed under the [MIT License](LICENSE).
Downloaded dependencies and the GCC compatibility patch remain subject to their
upstream licenses; see:
- [`f77_docker/THIRD_PARTY_NOTICES.md`](f77_docker/THIRD_PARTY_NOTICES.md)
- [`cpp_docker/THIRD_PARTY_NOTICES.md`](cpp_docker/THIRD_PARTY_NOTICES.md)
