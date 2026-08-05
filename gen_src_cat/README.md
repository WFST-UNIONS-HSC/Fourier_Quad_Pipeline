# External source catalog generators

This directory provides two ways to create the one-degree external source
catalog tiles consumed when the Fourier_Quad pipeline runs with `ext_cat = 1`:

- `process_extcat`: a C++17/MPI program that discovers existing raw text
  catalogs, projects their columns into the pipeline schema, and repartitions
  them by sky position.
- `query_y6gold_sync_mp_v2.py`: a Python program that queries the DES Y6 GOLD
  table through the NOIRLab Data Lab TAP service and writes the same tiles.

Both tools produce basenames such as
`des_y6_RA_299_300_Dec_m80_m79.dat` and the same 18-column commented-header
format. The C++ module is intentionally organized for a later merge into the
pipeline as `process_extcat`, alongside `process_init` and `process_main`; it is
not wired into `cpp_Standard` or `cpp_Lite` yet.

## C++ MPI catalog repartitioner

### Capabilities

The C++ implementation:

- reads any positive number of regular text catalogs from one directory;
- optionally scans subdirectories;
- accepts repeated, case-sensitive basename substring filters with OR
  semantics;
- detects whitespace, comma, or tab input independently for each file;
- recognizes canonical headers case-insensitively and reorders columns;
- accepts an explicit 18-column index projection for alternative or headerless
  schemas;
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
for the separately built tool. This ownership boundary matches the current
`process_init` and `process_main` execution model and must be preserved during
the later pipeline integration.

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
| `--columns LIST` | header names or first 18 fields | Exactly 18 comma-separated, one-based input indices in canonical output order |
| `--chunk-mib N` | `64` | Maximum nominal byte-range task size |
| `--malformed POLICY` | `fail` | `fail` or `skip` |
| `--existing POLICY` | `fail` | `fail` or `overwrite` generated tiles |

Boolean values accept `true`, `false`, `1`, `0`, `on`, or `off`. Named options
accept both `--name value` and `--name=value` forms.

### C++ input schema handling

With `--header auto`, each input file is inspected independently:

1. Leading blank lines and lines beginning with `#` are ignored as data.
2. A commented or uncommented header containing all 18 canonical field names
   is recognized case-insensitively.
3. Canonical fields may appear in any input order; they are projected into the
   fixed output order.
4. Without a recognized header, the first 18 fields are assumed to already be
   canonical.

Use `--columns` when an input uses different field names, contains extra fields,
or has no header. The list gives one-based input indices for the 18 canonical
output fields. For example, an input with an unused ID in column 1 followed by
the canonical values uses:

```bash
--header absent --columns 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19
```

All 18 projected tokens must be finite numbers because the downstream pipeline
reads them numerically. Quoted comma- or tab-delimited fields are accepted, but
embedded newlines are not. With `--malformed fail`, an invalid row stops the
collective job and no final tile set is published. With `--malformed skip`, the
row is counted and omitted.

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
merges shards in sorted input-file and original byte order, prepends one
canonical header, and publishes the complete tile set. Changing the MPI rank
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
catalogs and exercises one-, two-, and three-rank runs. It verifies canonical
header projection, explicit index projection, nested discovery, basename
filters, malformed-row skipping, RA/Dec boundaries, rank-count determinism,
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

## Shared output schema

Every final tile contains one commented header followed by 18 whitespace-
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

The pipeline reads `ra` and `dec` after skipping the number of leading fields
configured by `ext_cat_columns_before_ra`. Its default value is `4`, matching
the four flag fields above. External-catalog deblending expects the ten
magnitude/error fields followed by `dnf_z` and `dnf_zsigma` after `dec`, so the
generators always emit the complete canonical schema.

## Connect generated tiles to Fourier_Quad

Set `SOURCE_CAT` to the directory containing the generated `.dat` tiles, not to
an individual file:

- C++ Standard: `cpp_Standard/include/process_main/LensingConfig.hpp`
- C++ Lite: `cpp_Lite/include/process_main/LensingConfig.hpp`
- Fortran Standard: `f77/para.inc`
- Fortran Lite: `f77_Lite/para.inc`

For example, if the tiles are stored under `/data/catalogs/des_y6_chunks`, use
that directory as `SOURCE_CAT`, then rebuild the selected pipeline because this
setting is currently compiled into the executable.
