# Fourier_Quad_Pipeline

A weak-lensing shear measurement pipeline based on the **Fourier\_Quad** method,
provided in both **Fortran 77 (legacy)** and **C++ (current)** implementations,
together with ready-to-use **Docker** build environments and **HPC**
(Slurm/Apptainer) runner scripts.

> **中文文档**：请参阅 [README_CN.md](README_CN.md)

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

*The C++ version is recommended for modern environments.*

Each codebase has two variants:

- **Full** (`f77`, `cpp_Standard`): full feature set, including multi-scale /
  PCA PSF reconstruction (`PSFRecons` / `proc_psfreconsV2.f`).
- **Lite** (`f77_Lite`, `cpp_Lite`): frozen-branch simplified version. Eight
  compile-time switches are fixed to their production values and all dead-code
  branches are physically removed. See `cpp_Lite/REFACTOR_NOTES.md` for details.

---

## Quick Start

Two ways to run the pipeline. Choose a codebase (`f77` / `f77_Lite` or
`cpp_Standard` / `cpp_Lite`) and one of the modes below; see
[`F77_GUIDE.md`](F77_GUIDE.md) or [`CPP_GUIDE.md`](CPP_GUIDE.md) for full source,
configuration, Docker, and HPC-runner details.

First download the relevant source archives from the Releases.

### Mode A - Build from source

Manually install the toolchain, compile, and run the executable directly.

```bash
# Fortran 77: gfortran + MPICH + CFITSIO + LAPACK/BLAS
cd f77            # or f77_Lite
# edit para.inc (paths, PROCESS_stage, ...)
make
mpirun -np 4 ./Fourier_Quad_Pipe expo_list.list

# C++17: g++ + MPI + CFITSIO + FFTW3 + LAPACK/BLAS + Eigen3
cd cpp_Standard   # or cpp_Lite
# edit include/process_main/LensingConfig.hpp and include/ProcessConfig.hpp
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe --expo-list expo_list.list
```

Prerequisites and Makefile overrides: [`F77_GUIDE.md`](F77_GUIDE.md) /
[`CPP_GUIDE.md`](CPP_GUIDE.md). The full C++ parameter reference is
[`CPP_PIPELINE_PARAMETERS.md`](CPP_PIPELINE_PARAMETERS.md).

### Mode B - Use the runner with a prebuilt image

Suitable when you prefer not to configure the environment manually, or when
setting up the F77 environment is inconvenient.

Pull a published image from GHCR (or build it locally with the Docker
environment), then run locally or on an HPC cluster with the provided `runner/`
scripts.

```bash
# Pull a published image (or: cd f77_docker && docker compose build)
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/f77pipeline:latest
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/cpppipeline:latest

# HPC (Slurm + Apptainer): convert to a SIF and submit
bash f77_docker/runner/pull-sif.sh        # or cpp_docker/runner/pull-sif.sh
cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
# edit the .env paths for your cluster
sbatch f77_docker/runner/f77pipeline.slurm
```

Docker environment and HPC-runner details: [`F77_GUIDE.md`](F77_GUIDE.md) /
[`CPP_GUIDE.md`](CPP_GUIDE.md).

---

## Contributing

This project is an open-source implementation of Prof. Zhang Jun's series of
Fourier_Quad methods:

- Zhang, J. (2007). Measuring the cosmic shear in Fourier space: Measuring the cosmic shear in Fourier space. Monthly Notices of the Royal Astronomical Society, 383(1), 113–118. https://doi.org/10.1111/j.1365-2966.2007.12585.x
- Zhang, J., Luo, W., & Foucaud, S. (2015). Accurate shear measurement with faint sources. Journal of Cosmology and Astroparticle Physics, 2015(01), 024–024. https://doi.org/10.1088/1475-7516/2015/01/024
- Zhang, J., Zhang, P., & Luo, W. (2017). APPROACHING THE CRAMÉR–RAO BOUND IN WEAK LENSING WITH PDF SYMMETRIZATION. The Astrophysical Journal, 834(1), 8. https://doi.org/10.3847/1538-4357/834/1/8


## License

Repository-authored files are licensed under the [MIT License](LICENSE).
Downloaded dependencies and the GCC compatibility patch remain subject to their
upstream licenses; see:
- [`f77_docker/THIRD_PARTY_NOTICES.md`](f77_docker/THIRD_PARTY_NOTICES.md)
- [`cpp_docker/THIRD_PARTY_NOTICES.md`](cpp_docker/THIRD_PARTY_NOTICES.md)
