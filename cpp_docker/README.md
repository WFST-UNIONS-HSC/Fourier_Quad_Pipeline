# Fourier_Quad C++ container

This directory builds an x86_64 Linux toolchain image for `cpp_Standard` or
`cpp_Lite`. Source, catalogs, observation data, and outputs stay outside the
image on bind-mounted storage.

> 中文版：[README-CN.md](README-CN.md)

## Runtime

The Rocky Linux 8.10 image contains G++ 12.3.0, OpenMPI 4.1.8 with PMI2,
CFITSIO 4.6.4, FFTW 3.3.11, Eigen 3.4.0, LAPACK 3.11.0, and OpenBLAS 0.3.33.

Its portable HPC baseline is x86_64, Slurm `pmi2`, Apptainer/Singularity, a
shared filesystem, and routable TCP. Other architectures, PMIx-only sites,
schedulers, or vendor fabrics need separate validation.

## Build and verify

### Pull the GHCR image

```bash
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/cpppipeline:latest
```

### Download the source and build

Download the latest `cpp_docker.zip` from Releases.

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

## Local use

```bash
cp .env.example .env
# Set CPP_SOURCE_HOST and all host paths used by the selected phases.
docker compose run --rm FourierQuad-CPP
```

### Common `.env` parameters

Copy `.env.example`, then update the table for the host layout and selected
phases. `*_HOST` is a host path; `*_CONTAINER` is the absolute path visible to
the program inside the container.

| Parameter | Typical setting | Constraint |
|---|---|---|
| `IMAGE_NAME` | The image tag to run or build locally. | Must match the pulled image or `docker build -t` value. |
| `BUILD_JOBS` | Parallel jobs allowed while building the image. | Size for available CPU and memory. |
| `HOST_UID`, `HOST_GID` | Current host-user UID/GID. | Change when outputs must be directly writable by the host user. |
| `CPP_SOURCE_HOST` | The `cpp_Standard` or `cpp_Lite` source directory. | Mounted read/write at `/workspace/src_pipe` so build products persist. |
| `SCIENCE_ROOT_HOST/CONTAINER` | Science-image archive and its container path. | Needed only by `process_init`; use the container path in CLI arguments. |
| `DQ_ROOT_HOST/CONTAINER` | DQ-mask archive and its container path. | Required when DQ access is enabled; Lite always needs it. |
| `ASTROMETRY_CAT_HOST/CONTAINER` | Gaia catalog directory. | Container path must equal the compiled `ASTROMETRY_CAT`. |
| `SOURCE_CAT_HOST/CONTAINER` | Normalized External source catalog directory. | Container path must equal the effective `SOURCE_CAT` or `--extcat-output`. |
| `FLAT_PATH_HOST/CONTAINER` | Flat-calibration directory. | Needed only by a flat-enabled branch; container path must equal compiled `FLAT_PATH`. |
| `PROCESS_DATA_HOST/CONTAINER` | Writable processing directory. | Holds exposure lists, intermediates, and results; default container path is `/data/DataProcess`. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Set only when those phase paths need independent mounts. | Include `compose.optional.yaml` when used; otherwise use defaults below processing data. |

Inside the container:

```bash
make -C /workspace/src_pipe -j4
/workspace/src_pipe/Fourier_Quad_Pipe --help
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

Use container paths in program arguments. Core binds are source,
astrometry/source catalogs, flat calibration, and writable processing data.
Science/DQ archives and extcat/rearr/exposure-list/FD mounts are optional and
should be enabled only for phases that need them.

Catalog and calibration destinations used by compiled scientific branches
must match `config/LensingConfig.hpp`. `--extcat-output` can override the
external source-catalog tile path for one invocation.

For `process_astrocat`, expose the raw Gaia directory through a suitable
read-only bind and pass its container path with `--astrocat-input`. The
astrometry-catalog bind is read-only, so pass `--astrocat-output` to a writable
location, normally below `PROCESS_DATA_CONTAINER`. This option controls only
where the converter writes: it is not checked against and does not update the
compiled `LensingConfig::ASTROMETRY_CAT`. A later run that consumes the Type 2
tiles must bind them at the compiled astrometry path, set
`LensingConfig::AstroCatType = 2`, and rebuild separately.

## Slurm

Convert the same image to one SIF and follow the
[runner guide](runner/README.md). The supported launch boundary is:

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

See [SOURCES.md](SOURCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for dependency provenance.
