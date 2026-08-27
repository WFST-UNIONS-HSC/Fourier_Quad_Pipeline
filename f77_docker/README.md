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
