# cpp_Standard integrated workflow

`Fourier_Quad_Pipe` can repartition raw external source catalogs, prepare
compressed Science/DQ archives, run the numerical Fourier_Quad stages, and
spatially rearrange the generated `_all.cat` files in one MPI allocation. Each
phase can also run independently. Root `main.cpp` owns MPI initialization,
option parsing, phase ordering, and finalization. This full build retains PCA
`PSFRecons` support.

For the complete command-line, external-input, `ProcessConfig.hpp`, and
`LensingConfig.hpp` reference, see
[`CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Source layout

- `include/ProcessConfig.hpp`: workflow defaults and default phase switches.
- `include/process_extcat/`, `src/process_extcat/`: external-catalog schema,
  parsing, MPI byte-range partitioning, and deterministic tile publication.
- `include/process_init/`, `src/process_init/`: initializer wrapper plus the
  preserved `Initializer` and `FitsExtractor` modules.
- `include/process_main/`, `src/process_main/`: `LensingConfig`, all numerical
  modules, exposure-list loading, and the complete Stage 1–9 orchestration.
- `include/process_rearr/`, `src/process_rearr/`: self-contained `_all.cat`
  schema validation, full-sky weighted k-d partitioning, MPI redistribution,
  sorted subcatalog publication, and summary output.
- The `cpp_Standard` root contains only the executable entry point, build file,
  documentation, and phase implementation trees.

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

Run the focused external-catalog column reader test with:

```bash
make test-extcat-reader
```

Run the self-contained rearrangement unit and MPI integration tests with:

```bash
make test-rearr
```

## Defaults and option syntax

Edit `include/ProcessConfig.hpp` to set the normal dataset and default execution
mode. `RUN_PROCESS_EXTCAT`, `RUN_PROCESS_INIT`, and `RUN_PROCESS_REARR` default
to `false`; `RUN_PROCESS_MAIN` defaults to `true`. Every command-line option is
optional and overrides its configured default. Both `--name value` and
`--name=value` are accepted in any order.

The `EXTCAT_*` values configure raw-catalog discovery, the output directory,
parsing policies, MPI task size, optional ordered column selection, and raw
RA/Dec/ZP columns. With explicit selection disabled, every raw input field is
preserved in place. `process_main` reads only these three configured fields;
unselected catalog fields need not be numeric. `EXTCAT_TOTAL_COLUMNS` gives
`process_rearr` the pass-through external width; explicit projection instead
uses the projection-list length because that is the emitted external schema.
Set the primary catalog path with `SOURCE_CAT` in
`include/process_main/LensingConfig.hpp`. `EXTCAT_OUTPUT_DIRECTORY` is a
read-only reference to that value, so `process_extcat` writes where
`process_main` reads. `--extcat-output` overrides both for one invocation.

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
| `--run-extcat BOOL` | Enable or disable external-catalog repartitioning. |
| `--run-init BOOL` | Enable or disable archive initialization at runtime. |
| `--run-main BOOL` | Enable or disable the numerical pipeline at runtime. |
| `--run-rearr BOOL` | Enable or disable `_all.cat` spatial rearrangement; it follows `process_main` when both run. |
| `--extcat-input PATH` | Directory containing any number of raw text catalogs. |
| `--extcat-output PATH` | Destination tile directory and effective C++ `SOURCE_CAT`. |
| `--extcat-contains TEXT` | Repeatable case-sensitive basename token; repeats use OR. |
| `--extcat-recursive BOOL` | Enable or disable recursive input discovery. |
| `--extcat-delimiter MODE` | `auto`, `whitespace`, `comma`, or `tab`. |
| `--extcat-header MODE` | `auto`, `present`, or `absent`. |
| `--extcat-columns LIST` | One or more one-based raw indices; output fields follow this exact order and width. |
| `--extcat-ra-column N` | Raw one-based RA index; overrides `ra` header discovery. |
| `--extcat-dec-column N` | Raw one-based Dec index; overrides `dec` header discovery. |
| `--extcat-zp-column N` | Raw one-based photometric-redshift (`dnf_z`) index used by `process_main`. |
| `--extcat-chunk-mib N` | Positive MPI byte-range task size in MiB. |
| `--extcat-malformed POLICY` | `fail` or `skip` malformed rows. |
| `--extcat-existing POLICY` | `fail` or transactionally `overwrite` generated tiles. |
| `--science-root PATH` | Original read-only Science `.fits.fz` repository. |
| `--dq-root PATH` | Original read-only DQ `.fits.fz` repository. |
| `--output-root PATH` | Parent of the target directory and generated lists. |
| `--dataset TARGET:PREFIX` | One paired dataset; repeat the option for a batch. |
| `--target NAME` | Legacy single-dataset target; cannot be mixed with `--dataset`. |
| `--prefix TEXT` | Legacy single-dataset prefix; cannot be mixed with `--dataset`. |
| `--contains TEXT` | Accepted basename token; repeat for OR matching. |
| `--existing MODE` | `fail`, `resume`, or `overwrite`; default is `fail`. |
| `--f77-max-path N` | Maximum generated path length; `0` disables the check. |
| `--expo-list PATH` | Exposure list used in main/rearr-only mode. |
| `--help` | Print the effective command contract. |

Boolean values accept `true`/`false`, `1`/`0`, and `on`/`off`. One legacy
positional exposure-list path is retained as a compatibility alias, but new jobs
should use `--expo-list`. The first explicit `--dataset` replaces configured
`DATASETS`, and subsequent occurrences append. `--contains` follows the same
replacement/append rule for `CONTAINS`. Other duplicate scalar options use the
last value. The first explicit `--extcat-contains` similarly replaces
`EXTCAT_FILENAME_TOKENS`; later occurrences append. Dataset target names must be
unique within one invocation.

## Run modes

External-catalog-only execution accepts any number of matching input files and
does not require a configured dataset:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --extcat-input /data/raw_catalogs \
  --extcat-output /data/catalogs/des_y6_chunks \
  --extcat-contains .csv --extcat-contains y6_gold
```

`process_extcat` runs collectively once before the dataset loop. It writes the
same one-degree filenames as `gen_src_cat/query_y6gold_sync_mp_v2.py`. Its
default output preserves the complete raw schema; `--extcat-columns 5,3,4,1`,
for example, writes raw columns 5, 3, 4, and 1 as output columns 1–4. If it
fails, no later phase starts. Output passed onward to `process_main` must
include the raw columns configured by `EXTCAT_RA_COLUMN_ONE_BASED`,
`EXTCAT_DEC_COLUMN_ONE_BASED`, and `EXTCAT_ZP_COLUMN_ONE_BASED`; a
rearrangement-only run requires RA and Dec but does not require ZP. In
pass-through mode the required fields are read at their configured positions.
With explicit projection, the reader automatically maps each raw index to its
position in the ordered projection list.

Main-only local execution:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true \
  --expo-list /data/work/expo_g2019.list
```

Rearrangement-only local execution consumes existing per-exposure `_all.cat`
files referenced by the same exposure list:

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main false --run-rearr true \
  --expo-list /data/work/expo_g2019.list
```

All rearrangement-specific parameters are in
`include/process_rearr/ProcessRearrConfig.hpp`. With the default 18-field
external catalog, `ichi2=25` and the complete row width is calculated there as
`18 + 1 + 25 = 44`. The default 0.1-degree grid targets about 500,000 rows per
weighted k-d partition. Outputs are written below each dataset root in
`rearranged_catalog/` as `subcat_NNNNNN.cat` plus
`catalog_summary.txt`. Every data row must have the exact numeric width and
finite values; configured missing catalogs and malformed rows are skipped and
reported.

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

Chained local execution uses the same initializer options with downstream phase
switches enabled. After successful initialization, `process_main` and
`process_rearr` receive
the normalized absolute `output_root/expo_<target>.list` path returned by
`process_init`. That generated path overrides `--expo-list`, the legacy positional
argument, and every configured exposure-list default.

```bash
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init true --run-main true --run-rearr true \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work --dataset g2019:c4d_19 \
  --existing resume
```

Enable all four switches to build the external catalog first, then initialize,
process, and rearrange every configured dataset. The effective
`--extcat-output` path is also used by the numerical source extractor.

Datasets execute sequentially on the same MPI communicator and stop at the first
failure. In main/rearr-only batch mode, omit `--expo-list`; the driver derives one
`output_root/expo_<target>.list` path per dataset. A single external exposure list
is accepted only for a single downstream-only dataset. In chained batch mode, every
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
