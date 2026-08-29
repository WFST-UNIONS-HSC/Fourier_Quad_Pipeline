# Fourier_Quad Fortran on Slurm/Apptainer

This runner binds a writable Fortran source tree and data into an immutable
SIF, compiles once on the batch node, and launches one container command per
MPI rank.

> 中文版：[README-CN.md](README-CN.md)

## 1. Inspect the site

```bash
bash inspect-cluster-mpi.sh
```

Save the output and choose a launch mode only after checking MPI build options,
`srun --mpi=list`, Apptainer availability, and the compute-node fabric.

- `mpiexec` mode requires a host MPICH launcher compatible with the MPICH
  4.1.2 application in the SIF.
- `srun` mode requires a Slurm PMI interface explicitly validated with the SIF.
  The pilogin example uses `srun --mpi=pmi2`; host OpenMPI `mpirun` is never
  used to launch the MPICH application.

## 2. Configure shared paths

```bash
cp f77pipeline.env.example f77pipeline.env
```

Use `f77pipeline.pilogin-openmpi.env.example` only as a site-specific example,
not as a portable default.

### Common runner parameters

Every host path must be visible at the same location from all allocated nodes.
Paired container paths must match paths compiled into `para.inc` or the program's
positional argument.

| Parameter or location | Typical setting | Constraint |
|---|---|---|
| `OCI_IMAGE_URI`, `F77_SIF` | Remote F77 image and shared target/existing SIF path. | The SIF must be visible at the same path from every node. |
| `F77_SOURCE_HOST/CONTAINER` | `f77` or `f77_Lite`, normally mounted at `/workspace/f77`. | The host source must be writable so compilation products persist. |
| `ASTROMETRY_CAT_*`, `SOURCE_CAT_*`, `FLAT_PATH_*` | Gaia, External source catalog, and flat directories. | Container paths must match those compiled into `para.inc`. |
| `PROCESS_DATA_*`, `F77_EXPO_LIST_CONTAINER` | Shared writable processing directory and default container exposure-list path. | The list must resolve below a bind; an explicit positional argument overrides it. |
| `HPC_SHARED_SCRATCH_HOST` | Shared writable directory for the MPI smoke test. | Must be visible from all nodes and have sufficient space. |
| `APPTAINER_BIN`, `HPC_MODULES`, `SITE_ENV_SCRIPT` | Site Apptainer/Singularity command and module setup. | Keep `HPC_MODULES` as a Bash array. |
| `MPI_LAUNCH_MODE`, `MPI_LAUNCHER`, `SLURM_MPI_TYPE` | The generic template defaults to `mpiexec`; switch to `srun`/`pmi2` only after site validation. | The launcher must be compatible with MPICH 4.1.2 in the SIF; never use host OpenMPI `mpirun`. |
| `HPC_EXTRA_BINDS` | Comma-separated bind list for required extra host files/directories. | Preserve the template's string format; do not convert it to a Bash array. |
| `FI_PROVIDER`, `FI_PROVIDER_PATH` | Set only when site diagnostics require a specific libfabric provider. | Empty values retain the image MPICH/libfabric defaults. |
| `HPC_SCRUB_OPENMPI_ENV` | Normally keep `1`. | Prevents host OpenMPI transport variables from contaminating the MPICH container. |
| `F77_BUILD_JOBS`, `F77_MAKE_CLEAN`, `F77_EXECUTABLE` | Build parallelism, pre-build cleaning, and executable container path. | Do not let concurrent jobs clean or compile the same source copy. |
| `#SBATCH` in `f77pipeline.slurm` | Partition, nodes, tasks, tasks per node, CPUs, memory, time, and log paths. | Follow site policy and keep task counts consistent with the selected MPI mode. |

The env file is Bash. Keep `HPC_MODULES` as an indexed array. Container catalog
destinations must match the paths compiled into `para.inc`.

The mounted source Makefile must resolve libraries inside the SIF. This
repository's default library paths are site-specific; before submission,
adapt a private Makefile copy to `/opt/f77stack/lib` or otherwise ensure that a
bare `make -C "$F77_SOURCE_CONTAINER"` succeeds inside the image.

## 3. Acquire and verify the SIF

Set `OCI_IMAGE_URI` for the remote registry, then download the OCI image:

```bash
bash pull-sif.sh
```

`pull-sif.sh` refuses to overwrite an existing SIF. We then recommend running:

```bash
bash run-apptainer.sh --check
```

If login nodes cannot pull images, create the SIF on a permitted x86_64 Linux
host and transfer it to the configured shared path.

## 4. Validate MPI

Adjust the Slurm resource template and run:

```bash
sbatch mpi-smoke-test.slurm
```

For the pilogin PMI2 example, use the matching env file and
`mpi-smoke-test-pilogin-openmpi.slurm`. For a new site, we recommend running
single-rank, same-node multi-rank, and multi-node tests.

## 5. Run the pipeline

```bash
sbatch f77pipeline.slurm
```

With no script arguments, the executable receives
`F77_EXPO_LIST_CONTAINER`. An explicit positional list overrides it:

```bash
sbatch f77pipeline.slurm /data/DataProcess/another_expo.list
```

The job verifies binds, optionally cleans, compiles once, changes to the
processing directory, and launches all ranks. Do not let concurrent jobs clean
or compile the same source copy.

Resource directives, modules, and MPI mode are site inputs. Change them as a
coherent site configuration; do not assume that a validation from another
cluster proves compatibility with the local Slurm/PMI or interconnect.
