# Fourier_Quad Fortran container

This directory builds a reproducible x86_64 toolchain for `f77` or
`f77_Lite`. Source, catalogs, calibration files, processing data, and outputs
remain on bind-mounted host storage.

> 中文版：[README-CN.md](README-CN.md)

## Runtime

| Component | Version |
|---|---|
| Rocky Linux | 8.10 |
| GCC / GFortran | 4.8.5 |
| MPICH | 4.1.2 (`ch4:ofi`) |
| CFITSIO | 4.3.1 |
| LAPACK / reference BLAS | 3.8.0 |

Dependency sources and checksums are recorded in [SOURCES.md](SOURCES.md) and
`checksums.sha256`.

## Build and verify

### Pull the GHCR image

```bash
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/f77pipeline:latest
```

### Download the source and build

Download the latest `f77_docker.zip` from Releases.

```bash
docker build --platform linux/amd64 --build-arg BUILD_JOBS=4 \
  -t f77pipeline-dev:gnu4.8.5 .
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

## Local use

```bash
cp .env.example .env
# Set F77_SOURCE_HOST and the catalog/calibration/processing host paths.
docker compose run --rm FourierQuad-F77
```

### Common `.env` parameters

Copy `.env.example`, then update the table for the actual host layout.
`*_HOST` is a host path; `*_CONTAINER` is the absolute path visible to the
program inside the container.

| Parameter | Typical setting | Constraint |
|---|---|---|
| `IMAGE_NAME` | The F77 image tag to run or build locally. | Must match the pulled image or `docker build -t` value. |
| `BASE_IMAGE` | Base-image registry reference or pinned digest. | Change only when rebuilding the toolchain image. |
| `BUILD_JOBS` | Parallel jobs for building image dependencies. | Size for available CPU and memory. |
| `HOST_UID`, `HOST_GID` | Current host-user UID/GID. | Change when outputs must be directly writable by the host user. |
| `F77_SOURCE_HOST` | The `f77` or `f77_Lite` source directory. | Mounted read/write at `/workspace/f77` so build products persist. |
| `ASTROMETRY_CAT_HOST/CONTAINER` | Gaia catalog directory. | Container path must match the astrometry-catalog path in `para.inc`. |
| `SOURCE_CAT_HOST/CONTAINER` | External source catalog directory. | Container path must match the source-catalog path in `para.inc`. |
| `FLAT_PATH_HOST/CONTAINER` | Flat-calibration directory. | Needed only by a flat-enabled branch and must match `para.inc`. |
| `PROCESS_DATA_HOST/CONTAINER` | Writable processing directory. | Holds exposure lists, intermediates, and results; program arguments use container paths. |

Inside the container, compile the mounted source with the image libraries:

```bash
make -C /workspace/f77 clean
make -C /workspace/f77 \
  LAPACK_LIB_DIR=/opt/f77stack/lib \
  CFITSIO_LIB_DIR=/opt/f77stack/lib -j4
mpiexec -n 4 /workspace/f77/Fourier_Quad_Pipe \
  /data/DataProcess/expo_list.list
```

The repository Makefile has site-specific default library directories;
portable/container builds must override them as above or adapt a private
Makefile copy. Paths compiled into `para.inc` must be container paths matching
`ASTROMETRY_CAT_CONTAINER`, `SOURCE_CAT_CONTAINER`, and any active
`FLAT_PATH_CONTAINER` bind.

Rebuild the executable after changing source or any include file. Rebuild the
image only when the toolchain, dependencies, Dockerfile, or compatibility patch
changes.

## HPC

Docker Compose is for local use. On Slurm, convert the published image to a SIF
and follow the [runner guide](runner/README.md). The generic runner supports a
compatible host-MPICH `mpiexec` hybrid launch or an explicitly validated
Slurm-PMI launch. Never use host OpenMPI `mpirun` to launch the MPICH-linked
application.

Third-party licensing is summarized in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
