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

Inside the container:

```bash
make -C /workspace/src_pipe -j4
/workspace/src_pipe/Fourier_Quad_Pipe --help
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

Use container paths in program arguments. Core binds are source,
astrometry/source catalogs, flat calibration, and writable processing data.
Science/DQ archives and extcat/rearr/exposure-list/FD mounts are optional and
should be enabled only for phases that need them.

Catalog and calibration destinations used by compiled scientific branches
must match `config/LensingConfig.hpp`. `--extcat-output` can override the
external source-catalog tile path for one invocation.

## Slurm

Convert the same image to one SIF and follow the
[runner guide](runner/README.md). The supported launch boundary is:

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

See [SOURCES.md](SOURCES.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for dependency provenance.
