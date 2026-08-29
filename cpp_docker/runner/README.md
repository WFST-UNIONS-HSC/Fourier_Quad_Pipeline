# Fourier_Quad C++ on Slurm/Apptainer

This runner launches one precompiled SIF on x86_64 Slurm clusters that expose
`pmi2`. The application uses compiler, OpenMPI, and scientific libraries from
the SIF.

> 中文版：[README-CN.md](README-CN.md)

## Prerequisites and configuration

Confirm `pmi2` appears in `srun --mpi=list`, Apptainer/Singularity works on
compute nodes, and all paths are shared at identical locations. Then run:

```bash
bash inspect-cluster-mpi.sh
cp cpppipeline.env.example cpppipeline.env
```

### Common runner parameters

Host paths in `cpppipeline.env` must be visible at identical locations from all
allocated nodes. Paired container paths must match the program configuration or
CLI paths.

| Parameter or location | Typical setting | Constraint |
|---|---|---|
| `OCI_IMAGE_URI` / `CPP_DOCKER_ARCHIVE` | Set one according to OCI pulling or local Docker-archive conversion. | The URI names the remote image; the archive must be an x86_64 image visible to the build job. |
| `CPP_SIF`, `CPP_SIF_SHA256_EXPECTED` | Shared SIF path and optional expected SHA256. | The path must be identical on every node; fill the checksum after acquisition. |
| `CPP_SOURCE_HOST/CONTAINER` | `cpp_Standard` or `cpp_Lite`, normally mounted at `/workspace/src_pipe`. | The host source must be writable so compilation products persist. |
| `SCIENCE_ROOT_*`, `DQ_ROOT_*` | Science/DQ archives used by initialization. | Bind only when needed; Lite main processing requires per-chip DQ. |
| `ASTROMETRY_CAT_*`, `SOURCE_CAT_*`, `FLAT_PATH_*` | Gaia, External source catalog, and flat directories. | Container paths must match compiled paths or supported CLI overrides. |
| `PROCESS_DATA_*`, `CPP_EXPO_LIST_CONTAINER` | Shared writable processing directory and default container exposure-list path. | The exposure list must resolve below a bound container path. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Set only for phases needing independent directories. | Unset pairs are not bound; prefer defaults below `PROCESS_DATA`. |
| `HPC_SHARED_SCRATCH_HOST`, `APPTAINER_CACHE_DIR`, `APPTAINER_TMP_DIR` | Site shared writable scratch/cache/tmp paths. | Must have sufficient space and be compute-node accessible. |
| `APPTAINER_BIN`, `HPC_MODULES`, `SITE_ENV_SCRIPT` | Site command and module setup. | Keep `HPC_MODULES` as a Bash array and do not inject host MPI libraries. |
| `MPI_LAUNCH_MODE=srun`, `SLURM_MPI_TYPE=pmi2`, `SRUN_ARGS=()` | Preserve the runner launch mode; add site arguments through `SRUN_ARGS` if needed. | This runner requires `srun` + `pmi2`; do not flatten the array syntax. |
| `HPC_EXTRA_BINDS=()`, `HPC_PASSTHROUGH_ENV`, `HPC_CONTAINER_ENV=()` | Set only for required extra files or environment variables. | Keep Bash arrays and minimize host state entering the `--cleanenv` container. |
| `CPP_BUILD_JOBS`, `CPP_MAKE_CLEAN`, `CPP_EXECUTABLE` | Build parallelism, pre-build cleaning, and executable container path. | Do not let concurrent jobs clean or compile the same source copy. |
| `CPP_IMAGE_ID_EXPECTED` and expected compiler/MPI versions | Update together only when changing the image software stack. | Must match the versions actually present in the image. |
| `#SBATCH` in `cpppipeline.slurm` | Partition, nodes, tasks, tasks per node, CPUs, memory, time, and log paths. | Follow site policy and keep task counts consistent with the planned MPI ranks. |

The env file is Bash. Preserve its indexed arrays and keep
`MPI_LAUNCH_MODE=srun` with `SLURM_MPI_TYPE=pmi2`. Modules may expose Slurm or
Apptainer but must not inject host MPI libraries into the application.

## Acquire and validate

Use one acquisition path:

```bash
sbatch build-sif.slurm          # existing local Docker archive
bash pull-sif.sh                # download an OCI image from a remote registry
```

Both refuse to overwrite an existing SIF and create a SHA256 sidecar. We then
recommend running:

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

## Run

Pass C++ CLI options after the script name:

```bash
sbatch cpppipeline.slurm \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

With no arguments, the runner passes `CPP_EXPO_LIST_CONTAINER` as the legacy
exposure-list argument. The source bind is writable for compilation; do not
build the same copy concurrently. Adjust Slurm resource directives for the
site without changing the `srun --mpi=pmi2` launch boundary.
