# External source catalog generators

This directory provides two ways to create the one-degree external source
catalog tiles consumed when the Fourier_Quad pipeline runs with `ext_cat = 1`:

- `process_extcat`: a C++17/MPI program that discovers existing raw text
  catalogs, optionally selects columns in a user-specified order, and
  repartitions rows by sky position.
- `query_y6gold_sync_mp_v2.py`: a Python program that queries the DES Y6 GOLD
  table through the NOIRLab Data Lab TAP service and writes the same tiles.

Both tools produce basenames such as
`des_y6_RA_299_300_Dec_m80_m79.dat`. The Python downloader always writes its
DES Y6 GOLD 18-column query schema. The C++ repartitioner preserves the raw
schema by default and writes an arbitrary ordered subset only when explicit
column selection is enabled. It can reproduce the Python schema when the raw
catalog or selected list already has that order. The same reusable C++ module
is integrated into both `cpp_Standard` and `cpp_Lite` as the optional first
phase, before `process_init` and `process_main`. This directory retains the
independently buildable tool for catalog-only jobs and release packaging.

## C++ MPI catalog repartitioner

### Capabilities

The C++ implementation:

- reads any positive number of regular text catalogs from one directory;
- optionally scans subdirectories;
- accepts repeated, case-sensitive basename substring filters with OR
  semantics;
- detects whitespace, comma, or tab input independently for each file;
- preserves every input field in its original order when column selection is
  disabled;
- accepts any non-empty ordered list of one-based input indices when column
  selection is enabled, including repeated indices;
- locates raw `ra` and `dec` header fields independently from the output list,
  with explicit coordinate indices for headerless or nonstandard schemas;
- splits even one large input file into newline-aligned byte-range MPI tasks;
- writes one-degree sky tiles with deterministic row order independent of MPI
  rank count;
- supports fail-fast or skip behavior for malformed rows; and
- publishes output through isolated staging, with fail-on-existing or rollback-
  protected overwrite behavior.

The repartitioner does not deduplicate rows. If selected input catalogs overlap,
their repeated sources remain repeated in the output.

### Source layout and integration boundary

```text
gen_src_cat/
├── include/process_extcat/process_extcat.hpp
├── src/process_extcat/process_extcat.cpp
├── main.cpp
├── Makefile
└── tests/test_process_extcat.py
```

The reusable entry point is:

```cpp
int process_extcat(ProcessExtcat::Config config,
                   MPI_Comm communicator = MPI_COMM_WORLD);
```

`process_extcat` assumes MPI is already initialized and never calls
`MPI_Init` or `MPI_Finalize`. The standalone `main.cpp` owns those calls only
for the separately built tool. The copies under each C++ pipeline variant add a
`ProcessConfig::RuntimeOptions` adapter, while pipeline `main.cpp` retains MPI
lifecycle ownership.

### C++ build requirements

- C++17 compiler
- MPI implementation with a C++ compiler wrapper (`mpicxx`)
- GNU Make
- Python 3 only for the integration test suite

Local validation used GCC/mpicxx 15.2.0 and Open MPI 5.0.10 under Linux. The
production target is a Linux HPC cluster with a site-provided C++17 compiler
and MPI module. No CFITSIO, FFTW, Eigen, LAPACK, or other pipeline science
libraries are required by this standalone module.

Portable build:

```bash
cd gen_src_cat
make -j4
```

Override the MPI wrapper when a cluster uses a non-default command:

```bash
make CXX=/path/to/mpicxx -j4
```

Clean generated objects and the standalone executable with `make clean`.

### C++ quick start

The minimal command selects every regular file below the input directory:

```bash
mpirun -np 4 ./process_extcat \
  --input-dir /data/raw_catalogs \
  --output-dir /data/catalogs/des_y6_chunks
```

Select basenames containing `.csv` or `y6_gold`:

```bash
mpirun -np 4 ./process_extcat \
  --input-dir /data/raw_catalogs \
  --output-dir /data/catalogs/des_y6_chunks \
  --contains .csv \
  --contains y6_gold
```

Repeated `--contains` values use OR matching. In this example, a file is
selected when its basename contains either substring. Matching is
case-sensitive and applies to the basename, not the complete path.

On a Slurm cluster, use the site's supported launcher. A typical direct launch
is:

```bash
srun --ntasks=16 ./process_extcat \
  --input-dir /shared/raw_catalogs \
  --output-dir /shared/des_y6_chunks
```

### C++ command-line options

| Option | Default | Meaning |
|---|---|---|
| `--input-dir PATH` | required | Root containing raw catalogs |
| `--output-dir PATH` | required | Destination for final one-degree tiles |
| `--contains TEXT` | no filter | Repeatable basename substring; repeats use OR |
| `--recursive BOOL` | `true` | Scan below nested input directories |
| `--delimiter MODE` | `auto` | `auto`, `whitespace`, `comma`, or `tab` |
| `--header MODE` | `auto` | `auto`, `present`, or `absent` |
| `--columns LIST` | disabled | One or more comma-separated, positive one-based input indices; output fields follow this exact order |
| `--ra-column N` | named `ra`, otherwise `5` | Positive one-based raw RA index; overrides header-name discovery |
| `--dec-column N` | named `dec`, otherwise `6` | Positive one-based raw Dec index; overrides header-name discovery |
| `--chunk-mib N` | `64` | Maximum nominal byte-range task size |
| `--malformed POLICY` | `fail` | `fail` or `skip` |
| `--existing POLICY` | `fail` | `fail` or `overwrite` generated tiles |

Boolean values accept `true`, `false`, `1`, `0`, `on`, or `off`. Named options
accept both `--name value` and `--name=value` forms.

### C++ input schema handling

With `--header auto`, each input file is inspected independently. A unique pair
of case-insensitive `ra` and `dec` field names identifies a header immediately;
a leading nonnumeric record can also be recognized as a header while locating
the first valid coordinate row. Use `--header present` for a commented header
that uses nonstandard coordinate names, and `--header absent` for headerless
data.

Without `--columns`, every field is copied in its original order and each data
row must have the inspected input width. With `--columns`, the output width is
the list length and each item selects that one-based raw field. For example:

```bash
--columns 5,3,4,1
```

writes raw columns 5, 3, 4, and 1 as output columns 1, 2, 3, and 4. Indices may
repeat, in which case the raw field is repeated. Named input fields become the
output header in selected order; headerless fields use names such as
`column_5` and `column_3`.

Sky tiling always reads RA and Dec from the raw row, independently of the
output list. A recognized `ra`/`dec` header pair is preferred. Otherwise raw
columns 5 and 6 are the fallback. Use `--ra-column` and `--dec-column` together
for headerless or nonstandard layouts; either option enables explicit
coordinate indexing and therefore disables header-name discovery.

Only the raw RA and Dec tokens must be finite numbers. Other copied fields may
be strings. Quoted comma- or tab-delimited fields are accepted, but embedded
newlines are not. All selected input files must produce the same effective
output header; incompatible files fail before any final tile is published.
With `--malformed fail`, an invalid row stops the collective job and no final
tile set is published. With `--malformed skip`, the row is counted and omitted.

Valid right ascension is `[0, 360]` degrees; exactly 360 degrees wraps to zero.
Valid declination is `[-90, 90]` degrees; the exact north pole is placed in the
`p89_p90` tile. Other out-of-range coordinates are malformed.

The output directory may not equal or be nested below the input directory. This
prevents generated tiles or staging files from being rediscovered as raw input.

### C++ parallel and output behavior

Rank zero discovers and inspects sorted input paths. The data portion of every
file is divided into half-open byte ranges, with an effective chunk size small
enough to expose work to all ranks when the data volume permits. Each rank seeks
to assigned ranges and aligns to complete newline-delimited records.

Workers write collision-free per-task shards. After all tasks succeed, rank zero
merges shards in sorted input-file and original byte order, prepends the shared
projected header, and publishes the complete tile set. Changing the MPI rank
count therefore does not change output bytes.

The default `--existing fail` refuses to start when the output directory already
contains pipeline tile basenames. `--existing overwrite` moves the previous
generated tile set into private staging, publishes the new set, and restores the
old set if publication fails. Successful overwrite removes stale generated
tiles but preserves unrelated files in the output directory. If merge or
publication cannot complete, the error reports and retains the private staging
path so rollback artifacts are not destroyed before inspection.

### C++ tests

```bash
make test
```

The dependency-free integration test creates temporary CSV and whitespace
catalogs and exercises one-, two-, and three-rank runs. It verifies variable-
width pass-through, exact ordered `{5,3,4,1}` projection, headerless coordinate
indices, schema-mismatch rejection, nested discovery, basename filters,
malformed-row skipping, RA/Dec boundaries, rank-count determinism,
fail-on-malformed publication isolation, unsafe nested-output rejection,
fail-on-existing behavior, overwrite cleanup, and preservation of unrelated
files.

## Python DES Y6 GOLD TAP downloader

`query_y6gold_sync_mp_v2.py` queries `des_dr2.y6_gold` directly and is useful
when no local raw catalog is available.

### Python runtime requirements

- Linux or another environment with multiprocessing support
- Python 3.10 or newer (syntax checked locally with Python 3.12.13)
- NumPy (locally available as 2.5.0)
- PyVO 1.9.1 as the reference release; it installs the required Astropy and
  Requests dependencies (Astropy 8.0.0 was available during local syntax checks)
- Network access to `https://datalab.noirlab.edu/tap`

The script is interpreted and does not require compilation. On an HPC system,
load a site-provided Python 3 module before creating the virtual environment.
The exact module name is cluster-specific.

Create an isolated environment and install the two direct Python dependencies:

```bash
cd gen_src_cat
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pyvo
```

### Python configuration

The script currently uses constants rather than command-line options. Review
these values before starting a large download:

| Setting | Default | Purpose |
|---|---:|---|
| `OUT_DIR` | `des_y6_chunks` | Output directory, relative to the current working directory |
| `CONCURRENT_PROCESSES` | `4` | Number of simultaneous synchronous TAP queries |
| `TAP_URL` | NOIRLab Data Lab TAP endpoint | Catalog service URL |
| `MAX_ROWS` | `300000` | `TOP` limit applied to every query |
| `TARGET_RA_MIN`, `TARGET_RA_MAX` | `299`, `360` degrees | Half-open right-ascension range |
| `TARGET_DEC_MIN`, `TARGET_DEC_MAX` | `-80`, `20` degrees | Half-open declination range |
| `CHUNK_SIZE` | `1` degree | Width and height of each tile |

The sky-range and chunk-size constants are defined inside `main()`. Start with
low concurrency if the service is busy; reduce `CONCURRENT_PROCESSES` to 2 or 3
when requests are throttled or time out.

Keep integral-degree tile boundaries unless the filename formatter and the
pipeline lookup logic are updated together.

### Run the Python downloader

Run from this directory so that the default output is created at
`gen_src_cat/des_y6_chunks/`:

```bash
python query_y6gold_sync_mp_v2.py
```

Each worker performs a blocking TAP query. Successful non-empty results are
written in Astropy's `ascii.commented_header` format. Existing tiles containing
fewer than `MAX_ROWS` data rows are skipped, so rerunning the command resumes
the grid after failures.

### Python row-limit and retry behavior

A result with exactly `MAX_ROWS` rows may have been truncated by the ADQL `TOP`
clause. The script prints a warning in that case. On a later run it deletes the
tile and repeats the same query; it does **not** automatically subdivide the
region. Treat such tiles as incomplete until the query is manually partitioned
with a compatible naming and lookup strategy or the service can return the full
result.

Network errors are reported per tile and do not stop the entire process pool.
Rerun the script to retry tiles that were not written.

## DES Y6 GOLD schema and pipeline compatibility

The Python downloader writes one commented header followed by 18 whitespace-
delimited columns in this order:

| Column | Field |
|---:|---|
| 1 | `flags_footprint` |
| 2 | `flags_foreground` |
| 3 | `flags_gold` |
| 4 | `ext_mash` |
| 5 | `ra` |
| 6 | `dec` |
| 7 | `bdf_mag_g` |
| 8 | `bdf_mag_err_g` |
| 9 | `bdf_mag_r` |
| 10 | `bdf_mag_err_r` |
| 11 | `bdf_mag_i` |
| 12 | `bdf_mag_err_i` |
| 13 | `bdf_mag_z` |
| 14 | `bdf_mag_err_z` |
| 15 | `bdf_mag_y` |
| 16 | `bdf_mag_err_y` |
| 17 | `dnf_z` |
| 18 | `dnf_zsigma` |

The C++ repartitioner emits this schema only when its pass-through input or
explicit selection has the same fields in the same order. Otherwise it emits
the effective input width or selected-list width.

`process_main` reads only the one-based raw positions configured by
`EXTCAT_RA_COLUMN_ONE_BASED`, `EXTCAT_DEC_COLUMN_ONE_BASED`, and
`EXTCAT_ZP_COLUMN_ONE_BASED`; their defaults are `5`, `6`, and `17` for the
schema above. Other fields are not converted and may contain arbitrary strings.
When explicit projection is enabled, all three raw indices must appear in the
selected list. The reader automatically maps them to their generated output
positions, so no fixed width or fixed magnitude/error layout is required.

## Connect generated tiles to Fourier_Quad

For the C++ programs, configure `SOURCE_CAT` in the selected
`include/process_main/LensingConfig.hpp`, or pass `--extcat-output`.
`EXTCAT_OUTPUT_DIRECTORY` follows that primary path, and a command-line
override updates the effective `SOURCE_CAT` before either phase starts. Enable
`RUN_PROCESS_EXTCAT` or pass `--run-extcat true` to generate the tiles as the
first phase. The Fortran programs still configure `SOURCE_CAT` directly:

- Fortran Standard: `f77/para.inc`
- Fortran Lite: `f77_Lite/para.inc`

The integrated C++ command-line interface uses the same policies as the
standalone tool, prefixed with `--extcat-`; for example, standalone `--contains`
becomes `--extcat-contains`. See the selected C++ README for the complete
three-phase command contract.
