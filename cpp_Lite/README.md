# cpp_Lite integrated workflow

`Fourier_Quad_Pipe` can prepare compressed Science/DQ archives, run the numerical
Fourier_Quad stages, or run both phases in one MPI allocation. The initializer
and numerical pipeline retain separate source/header trees; root `main.cpp` owns
MPI initialization, option parsing, phase ordering, and finalization. This Lite
build preserves its frozen scientific branches and intentionally omits PCA
`PSFRecons` support.

## Source layout

- `include/ProcessConfig.hpp`: workflow defaults and default phase switches.
- `include/process_init/`, `src/process_init/`: initializer wrapper plus the
  preserved `Initializer` and `FitsExtractor` modules.
- `include/process_main/`, `src/process_main/`: `LensingConfig`, all numerical
  modules, exposure-list loading, and the complete Stage 1–9 orchestration.
- The `cpp_Lite` root contains only the executable entry point, build file,
  documentation, and the two implementation trees.

## Compiler and libraries

Local verification used GCC/G++ 15.2.0, Open MPI 5.0.10, CFITSIO 4.6.3, FFTW
3.3.10, Eigen 3.4.0, and OpenBLAS/LAPACK 0.3.33 on Linux in WSL2.

The existing portable HPC target is GCC/G++ 12.3.0, Open MPI 4.1.8, CFITSIO
4.6.4, FFTW 3.3.11, Eigen 3.4.0, and LAPACK 3.11.0. A cluster may use equivalent
site modules as long as one matching MPI C++ wrapper compiles and launches the
program.

## Build

When all headers and libraries are visible through the compiler's normal search
paths:

```bash
make -j4
```

For one consolidated scientific-stack prefix containing headers, Eigen under
`include/eigen3`, and libraries under `lib`:

```bash
make STACK_PREFIX="${STACK_PREFIX}" -j4
```

When Eigen and the linked libraries use different prefixes, pass portable local
overrides explicitly:

```bash
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" -j4
```

No Windows-native compiler or wrapper is required. Cluster builds use the same
Makefile after loading the site's compiler, MPI, and scientific-library modules.

## Defaults and option syntax

Edit `include/ProcessConfig.hpp` to set the normal dataset and default execution
mode. `RUN_PROCESS_INIT` defaults to `false`; `RUN_PROCESS_MAIN` defaults to
`true`. Every command-line option is optional and overrides its configured
default. Both `--name value` and `--name=value` are accepted in any order.

`DATASETS` stores paired target/prefix values, and `CONTAINS` stores the archive
basename tokens accepted with OR semantics. For example:

```cpp
inline const std::vector<DatasetSpec> DATASETS = {
    {"g2013", "c4d_13"},
    {"g2014", "c4d_14"},
    {"g2019", "c4d_19"},
};
inline const std::vector<std::string> CONTAINS = {"v1", "v2"};
```

| Option | Purpose |
|:---|:---|
| `--run-init BOOL` | Enable or disable archive initialization at runtime. |
| `--run-main BOOL` | Enable or disable the numerical pipeline at runtime. |
| `--science-root PATH` | Original read-only Science `.fits.fz` repository. |
| `--dq-root PATH` | Original read-only DQ `.fits.fz` repository. |
| `--output-root PATH` | Parent of the target directory and generated lists. |
| `--dataset TARGET:PREFIX` | One paired dataset; repeat the option for a batch. |
| `--target NAME` | Legacy single-dataset target; cannot be mixed with `--dataset`. |
| `--prefix TEXT` | Legacy single-dataset prefix; cannot be mixed with `--dataset`. |
| `--contains TEXT` | Accepted basename token; repeat for OR matching. |
| `--existing MODE` | `fail`, `resume`, or `overwrite`; default is `fail`. |
| `--f77-max-path N` | Maximum generated path length; `0` disables the check. |
| `--expo-list PATH` | Exposure list used in main-only mode. |
| `--help` | Print the effective command contract. |

Boolean values accept `true`/`false`, `1`/`0`, and `on`/`off`. One legacy
positional exposure-list path is retained as a compatibility alias, but new jobs
should use `--expo-list`. The first explicit `--dataset` replaces configured
`DATASETS`, and subsequent occurrences append. `--contains` follows the same
replacement/append rule for `CONTAINS`. Other duplicate scalar options use the
last value. Dataset target names must be unique within one invocation.

## Run modes

Main-only local execution:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true \
  --expo-list /data/work/expo_g2019.list
```

Initializer-only local execution:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work --dataset g2019:c4d_19
```

Batch initialization pairs every target with its own prefix. Multiple contains
tokens select an archive when any token occurs in its basename:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main false \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work \
  --dataset g2013:c4d_13 --dataset g2014:c4d_14 \
  --dataset g2019:c4d_19 \
  --contains v1 --contains v2
```

Chained local execution uses the same initializer options with both phase
switches enabled. After successful initialization, `process_main` always receives
the normalized absolute `output_root/expo_<target>.list` path returned by
`process_init`. That generated path overrides `--expo-list`, the legacy positional
argument, and every configured exposure-list default.

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main true \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work --dataset g2019:c4d_19 \
  --existing resume
```

Datasets execute sequentially on the same MPI communicator and stop at the first
failure. In main-only batch mode, omit `--expo-list`; the driver derives one
`output_root/expo_<target>.list` path per dataset. A single external exposure list
is accepted only for a single main-only dataset. In chained batch mode, every
initializer-generated absolute list overrides external exposure-list input for
its corresponding dataset.

On a Slurm cluster, use the same executable arguments with the site launcher,
for example `srun -n 40 ./Fourier_Quad_Pipe ...`. The initializer failure status
is collective; the numerical phase is never entered after initialization fails.

## Initializer output contract

For each `--output-root OUTPUT --dataset TARGET:PREFIX`, initialization creates the ten
pipeline subdirectories below `OUTPUT/TARGET`, writes uncompressed Science/DQ
chip FITS files, and publishes:

- `OUTPUT/expo_TARGET.list`
- `OUTPUT/fits_TARGET.list`
- `OUTPUT/init_TARGET_manifest.json`

Manifest schema version 2 records all active basename filters in the
`filename_tokens` array.

Science chips are numbered by two-dimensional HDU occurrence. DQ chips are
numbered by `CCDNUM`. Source `.fits.fz` archives are read in place and are never
copied or removed.
