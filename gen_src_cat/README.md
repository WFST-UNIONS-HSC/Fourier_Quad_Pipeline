# External source catalogs

This directory provides two user-facing tools:

- `process_extcat`: C++17/MPI repartitioner for existing text catalogs;
- `query_y6gold_sync_mp_v2.py`: multiprocessing TAP downloader for DES Y6
  GOLD when no local raw catalog is available.

The C++ repartitioner produces canonical one-degree pipeline tiles using
`PathConfig::SOURCE_CAT_TILE_PREFIX` from `cpp_Standard/config/pathconfig.hpp`.
Its default basename is `extern_RA_299_300_Dec_m80_m79.dat`; the configured
prefix excludes the fixed `RA_` token. The downloader writes raw source files,
which should be passed through `process_extcat` before pipeline use.

## MPI repartitioner

The standalone program requires only an MPI C++17 wrapper and GNU Make.

```bash
make -j4
mpirun -np 4 ./process_extcat \
  --input-dir /data/raw_catalogs \
  --output-dir /data/catalogs/des_y6_tiles
```

Current Make targets are `all` and `clean`.

### Options

| Option | Default | Meaning |
|---|---|---|
| `--input-dir PATH` | required | Raw catalog root. |
| `--output-dir PATH` | required | Final tile directory. |
| `--contains TEXT` | no filter | Repeatable basename substring; OR semantics. |
| `--recursive BOOL` | `true` | Scan subdirectories. |
| `--delimiter MODE` | `auto` | `auto`, `whitespace`, `comma`, or `tab`. |
| `--header MODE` | `auto` | `auto`, `present`, or `absent`. |
| `--columns LIST` | pass-through | Ordered, comma-separated one-based projection. |
| `--ra-column N` | named `ra`, else `5` | Explicit raw one-based RA field. |
| `--dec-column N` | named `dec`, else `6` | Explicit raw one-based Dec field. |
| `--chunk-mib N` | `64` | Nominal byte-range task size. |
| `--malformed POLICY` | `fail` | `fail` or `skip`. |
| `--existing POLICY` | `fail` | `fail` or `overwrite`. |

Use `./process_extcat --help` for the current interface. Options accept
`--name value` and `--name=value`; booleans accept `true/false`, `1/0`, and
`on/off`.

Without `--columns`, every raw field is preserved. With projection, output
width and order follow the list. Tiling always uses raw RA/Dec, independently
of the selected output columns. All input files must resolve to one compatible
output schema.

The program reads newline-aligned MPI byte ranges and publishes tiles in
deterministic input order. It does not deduplicate overlapping catalogs. The
output directory cannot equal or be nested below the input. `fail` policies
leave no partially published final tile set; `overwrite` stages the old set so
it can be restored if publication fails.

## DES Y6 GOLD TAP downloader

The Python script uses NumPy and PyVO:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pyvo
python query_y6gold_sync_mp_v2.py
```

It has no CLI. Review these constants in the script before a large download:

| Setting | Default |
|---|---:|
| `OUT_DIR` | `des_y6_chunks` |
| `CONCURRENT_PROCESSES` | `4` |
| `MAX_ROWS` | `300000` per query |
| RA range | `[299, 360)` degrees |
| Dec range | `[-80, 20)` degrees |
| `CHUNK_SIZE` | `1` degree |

The script queries the NOIRLab Data Lab TAP endpoint synchronously. Existing
tiles below `MAX_ROWS` are skipped, so rerunning resumes missing work. A result
equal to `MAX_ROWS` may be truncated: the script warns but does not subdivide
the tile automatically. Network failures are reported per tile and retried by
rerunning the script.

The downloader writes this 18-field DES schema:

```text
flags_footprint flags_foreground flags_gold ext_mash ra dec
bdf_mag_g bdf_mag_err_g bdf_mag_r bdf_mag_err_r
bdf_mag_i bdf_mag_err_i bdf_mag_z bdf_mag_err_z
bdf_mag_y bdf_mag_err_y dnf_z dnf_zsigma
```

## Use tiles in the pipeline

For C++, either run the integrated `process_extcat` phase with
`--run-extcat true` or use the standalone tool and pass its destination as
`--extcat-output`. The integrated option names add the `--extcat-` prefix to
standalone policies. The standalone executable uses the Standard
`pathconfig.hpp` by default. Build with `make PIPELINE_VARIANT=cpp_Lite` when it
must follow the Lite configuration instead; the selected producer and
SOURCE_CAT consumer then share exactly one `SOURCE_CAT_TILE_PREFIX`. Rebuild
after changing that value.

For Fortran, set `SOURCE_CAT` in `f77/para.inc` or `f77_Lite/para.inc` and
rebuild. The default pipeline schema expects RA, Dec, and photo-z at one-based
positions 5, 6, and 17.
