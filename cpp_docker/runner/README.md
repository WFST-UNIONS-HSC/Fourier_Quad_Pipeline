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

Set the image/SIF, source, processing, catalog, calibration, cache, and tmp
paths. Optional Science/DQ/extcat/rearr/exposure-list/FD binds are needed only
for the selected phases.

The env file is Bash. Preserve its indexed arrays and keep
`MPI_LAUNCH_MODE=srun` with `SLURM_MPI_TYPE=pmi2`. Modules may expose Slurm or
Apptainer but must not inject host MPI libraries into the application.

## Acquire and validate

Use one acquisition path:

```bash
sbatch build-sif.slurm          # reviewed Docker archive
bash pull-sif.sh                # digest-pinned OCI image
```

Both refuse to overwrite an existing SIF and create a SHA256 sidecar. Then:

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
