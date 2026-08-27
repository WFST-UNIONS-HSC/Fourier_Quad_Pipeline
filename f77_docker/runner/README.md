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
not as a portable default. Set the SIF, source, catalogs, calibration,
processing, exposure-list, scratch, and optional cache/runtime paths. Every
host path must be visible at the same location from all allocated nodes.

The env file is Bash. Keep `HPC_MODULES` as an indexed array. Container catalog
destinations must match the paths compiled into `para.inc`.

The mounted source Makefile must resolve libraries inside the SIF. This
repository's default library paths are site-specific; before submission,
adapt a private Makefile copy to `/opt/f77stack/lib` or otherwise ensure that a
bare `make -C "$F77_SOURCE_CONTAINER"` succeeds inside the image.

## 3. Acquire and verify the SIF

Set a digest-pinned `OCI_IMAGE_URI` when possible, then:

```bash
bash pull-sif.sh
bash run-apptainer.sh --check
```

`pull-sif.sh` refuses to overwrite an existing SIF. If login nodes cannot pull
images, create the SIF on a permitted x86_64 Linux host and transfer it to the
configured shared path.

## 4. Validate MPI

Adjust the Slurm resource template and run:

```bash
sbatch mpi-smoke-test.slurm
```

For the pilogin PMI2 example, use the matching env file and
`mpi-smoke-test-pilogin-openmpi.slurm`. A new site is ready only after
single-rank, same-node multi-rank, and multi-node tests succeed.

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
