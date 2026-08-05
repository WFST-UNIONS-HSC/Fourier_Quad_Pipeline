# C++ Pipeline External Inputs and Parameter Reference

This is the authoritative reference for inputs accepted by the current
`cpp_Standard` and `cpp_Lite` executables. It covers:

- runtime command-line options parsed by `main.cpp`;
- workflow defaults in `include/ProcessConfig.hpp`;
- scientific compile-time parameters in
  `include/process_main/LensingConfig.hpp`;
- external files and directories read by `process_extcat`, `process_init`, and
  `process_main`.

It does not describe the Fortran parameter files or Makefile/build variables.
MPI process count and launcher options belong to `mpirun`, `mpiexec`, or the
site scheduler and are not parsed by `Fourier_Quad_Pipe`.

## Configuration layers and precedence

The executable has two configuration layers:

1. `ProcessConfig.hpp` provides workflow defaults. A matching command-line
   option overrides the default for one invocation; rebuilding is not required.
2. `LensingConfig.hpp` provides scientific and numerical constants. Changes to
   this file require rebuilding the selected C++ variant.

`LensingConfig::SOURCE_CAT` is the primary compiled default for the external
catalog directory. `ProcessConfig::EXTCAT_OUTPUT_DIRECTORY` references it, and
`--extcat-output` updates the effective `SOURCE_CAT` for one invocation.

Command-line options accept both `--name value` and `--name=value`. Boolean
values accept `true`, `false`, `1`, `0`, `on`, and `off`. Duplicate scalar
options use the last value. The first explicit `--dataset`, `--contains`, or
`--extcat-contains` clears its configured list; later occurrences append.

`process_extcat` runs once, followed by `process_init` and `process_main` for
each dataset. At least one phase must be enabled. A failure in an earlier phase
prevents later phases from starting.

## Integrated C++ entry-point parameters

The root driver calls the three pipeline functions in the following order. It
owns `MPI_Init` and `MPI_Finalize`; none of these functions may initialize or
finalize MPI independently.

| Function | Input / output parameter | Meaning |
|:---|:---|:---|
| `process_extcat` | `const ProcessConfig::RuntimeOptions& options` | Integrated extcat settings described in the `RuntimeOptions` table below. |
| `process_extcat` | `MPI_Comm communicator` | Collective communicator; defaults to `MPI_COMM_WORLD`. Every participating rank must call the function with consistent settings. |
| `process_init` | `const ProcessConfig::RuntimeOptions& options` | Initializer roots, filters, existing-output policy, and path limit. |
| `process_init` | `const ProcessConfig::DatasetSpec& dataset` | One target/prefix pair from the effective dataset list. |
| `process_init` | `std::string& generated_expo_list` | Output parameter receiving the normalized absolute exposure-list path on success. |
| `process_main` | `const std::string& exposure_list` | Exposure-list file used to load and broadcast numerical pipeline jobs. |

The reusable lower-level extcat overload accepts a `ProcessExtcat::Config` plus
an optional communicator. Its fields map to integrated settings as follows:

| `ProcessExtcat::Config` field | Default | Integrated source |
|:---|:---|:---|
| `input_directory` | Empty | `RuntimeOptions::extcat_input_directory` |
| `output_directory` | Empty | `RuntimeOptions::extcat_output_directory` |
| `filename_tokens` | Empty | `RuntimeOptions::extcat_filename_tokens` |
| `recursive` | `true` | `RuntimeOptions::extcat_recursive` |
| `delimiter` | `Delimiter::Auto` | `RuntimeOptions::extcat_delimiter` |
| `header_mode` | `HeaderMode::Auto` | `RuntimeOptions::extcat_header_mode` |
| `malformed_policy` | `MalformedPolicy::Fail` | `RuntimeOptions::extcat_malformed_policy` |
| `existing_policy` | `ExistingPolicy::Fail` | `RuntimeOptions::extcat_existing_policy` |
| `chunk_bytes` | `64 * 1024 * 1024` | `RuntimeOptions::extcat_chunk_mib` converted to bytes |
| `use_explicit_columns` | `false` | `RuntimeOptions::extcat_use_explicit_columns` |
| `input_columns` | Zero-based `{0, ..., 17}` | Ordered list used only when `use_explicit_columns=true`; one-based `RuntimeOptions::extcat_input_columns_one_based` is converted to zero-based indices |
| `use_explicit_coordinate_columns` | `false` | `RuntimeOptions::extcat_use_explicit_coordinate_columns` |
| `ra_column` | Zero-based `4` | One-based `RuntimeOptions::extcat_ra_column_one_based` converted to zero-based; named `ra` overrides it unless explicit coordinate mode is enabled |
| `dec_column` | Zero-based `5` | One-based `RuntimeOptions::extcat_dec_column_one_based` converted to zero-based; named `dec` overrides it unless explicit coordinate mode is enabled |

`ProcessExtcat::kCanonicalColumnCount` and `canonicalColumnNames()` retain the
DES Y6 GOLD 18-field reference schema for API compatibility. They do not set
the generated width.

## External input artifacts

| Input | Configured by | Required when | Format and purpose |
|:---|:---|:---|:---|
| Raw external catalogs | `EXTCAT_INPUT_DIRECTORY` / `--extcat-input` | `process_extcat` is enabled | Any number of regular text files. Discovery may be recursive and may use case-sensitive basename substring filters. Supported delimiters are whitespace, comma, and tab. |
| External catalog tile directory | `SOURCE_CAT` / `EXTCAT_OUTPUT_DIRECTORY` / `--extcat-output` | `process_extcat` or `process_main` is enabled | `SOURCE_CAT` is the primary configured path. `EXTCAT_OUTPUT_DIRECTORY` follows it, and a CLI override updates the effective path before `process_extcat` or `process_main` starts. |
| Science archive repository | `SCIENCE_ROOT` / `--science-root` | `process_init` is enabled | Read-only multi-HDU Science FITS/FZ archives selected by dataset prefix and filename tokens. |
| DQ archive repository | `DQ_ROOT` / `--dq-root` | `process_init` is enabled | Read-only multi-HDU DQ FITS/FZ archives paired with the selected Science archives. |
| Pipeline output root | `OUTPUT_ROOT` / `--output-root` | `process_init`, or derived-list main-only execution | Parent of dataset directories and `expo_<target>.list`. It also contains the extracted chip images and generated intermediate products. |
| Exposure list | `EXPO_LIST`, `--expo-list`, positional `LEGACY_EXPO_LIST`, or the initializer result | `process_main` is enabled | Text file whose non-empty lines identify per-exposure chip-list files. An initializer-generated list takes precedence in chained execution. |
| Gaia tile catalog | `ASTROMETRY_CAT` | Gaia astrometry is selected | Directory of the Gaia reference tiles expected by `generateGaiaFileName`. |
| Legacy flat/mask files | `FLAT_PATH` | Standard build with `include_FLAT=1` or `include_Mask=1/3` | Per-chip FITS files used for super-flat multiplication and/or the legacy mask branch. |
| External PSF image | `PSF_PATH` | Standard build with `ext_PSF=1` | A FITS image named `PSF.fits` below `PSF_PATH`. |
| Per-chip DQ masks | Produced by `process_init` below `OUTPUT_ROOT/<target>/dqmask` | `include_Mask=2/3`; Lite is fixed to `2` | Uncompressed per-chip FITS masks consumed during preprocessing. |

All paths are resolved by the C++ process. The external-catalog output directory
must not equal, or be nested below, its input directory.

## Runtime command line

### Phase selection and external-catalog repartitioning

| Option | Config default | Accepted values | Function and constraints |
|:---|:---|:---|:---|
| `--run-extcat BOOL` | `RUN_PROCESS_EXTCAT=false` | Boolean | Run `process_extcat` as the first phase. Requires a non-empty input and output directory. |
| `--run-init BOOL` | `RUN_PROCESS_INIT=false` | Boolean | Run archive discovery/extraction and publish pipeline lists. |
| `--run-main BOOL` | `RUN_PROCESS_MAIN=true` | Boolean | Run the numerical Stage 1–9 pipeline. |
| `--extcat-input PATH` | `EXTCAT_INPUT_DIRECTORY=""` | Existing directory | Root containing raw catalog files. Required when `--run-extcat=true`. |
| `--extcat-output PATH` | `SOURCE_CAT` through `EXTCAT_OUTPUT_DIRECTORY` | Directory path | Override both the generated-tile destination and effective `SOURCE_CAT` for this invocation. Required when extcat or main runs. |
| `--extcat-contains TEXT` | Empty list | Repeatable non-empty text | Case-sensitive substring match against the basename only. A file matches if any token matches; an empty list accepts all regular files. |
| `--extcat-recursive BOOL` | `true` | Boolean | Recurse into subdirectories when discovering raw catalogs. |
| `--extcat-delimiter MODE` | `auto` | `auto`, `whitespace`, `comma`, `tab` | Select the raw table delimiter. `auto` chooses comma when present, otherwise tab when present, otherwise whitespace. |
| `--extcat-header MODE` | `auto` | `auto`, `present`, `absent` | Control header handling. `auto` recognizes unique case-insensitive `ra`/`dec` names and can classify a leading record as a header while finding the first coordinate row; `present` requires a header; `absent` treats records as headerless. Leading blank and unrecognized `#` comment lines are skipped. |
| `--extcat-columns LIST` | Disabled; stored list 1–18 | One or more positive, comma-separated, one-based indices | Enable explicit projection. Output column 1 uses the first listed raw field, output column 2 the second, and so on. List length sets output width; repeated indices are allowed. Without this option every raw field is preserved in place. |
| `--extcat-ra-column N` | Named `ra`, otherwise configured fallback | Positive one-based index, distinct from Dec | Select the raw RA field used only for sky tiling and enable explicit coordinate indexing. This overrides header-name discovery. |
| `--extcat-dec-column N` | Named `dec`, otherwise configured fallback | Positive one-based index, distinct from RA | Select the raw Dec field used only for sky tiling and enable explicit coordinate indexing. This overrides header-name discovery. |
| `--extcat-chunk-mib N` | `64` | Positive integer | Approximate newline-aligned MPI byte-range task size in MiB. It controls task granularity, not the final tile size. |
| `--extcat-malformed POLICY` | `fail` | `fail`, `skip` | Stop collectively on the first malformed data row, or skip malformed rows and report the count. |
| `--extcat-existing POLICY` | `fail` | `fail`, `overwrite` | Reject existing generated tiles, or transactionally replace the complete generated tile set. |

For example, `--extcat-columns 5,3,4,1` writes raw fields 5, 3, 4, and
1 as output fields 1–4. Only the raw RA and Dec fields must be finite numeric
values; projected payloads may be strings. All selected input files must have
the same effective projected header, or the job fails before publication.

The Python downloader's 18-column DES Y6 GOLD reference order is:

| Position | Name | Position | Name |
|---:|:---|---:|:---|
| 1 | `flags_footprint` | 10 | `bdf_mag_err_r` |
| 2 | `flags_foreground` | 11 | `bdf_mag_i` |
| 3 | `flags_gold` | 12 | `bdf_mag_err_i` |
| 4 | `ext_mash` | 13 | `bdf_mag_z` |
| 5 | `ra` | 14 | `bdf_mag_err_z` |
| 6 | `dec` | 15 | `bdf_mag_y` |
| 7 | `bdf_mag_g` | 16 | `bdf_mag_err_y` |
| 8 | `bdf_mag_err_g` | 17 | `dnf_z` |
| 9 | `bdf_mag_r` | 18 | `dnf_zsigma` |

Rows are assigned to integral-degree RA/Dec tiles using raw coordinate fields
that are independent of the output projection. The filename convention matches
`gen_src_cat/query_y6gold_sync_mp_v2.py`; the schema matches only when the raw
input or explicit list has the 18 fields above in that order.

### Initializer and main-pipeline inputs

| Option | Config default | Accepted values | Function and constraints |
|:---|:---|:---|:---|
| `--science-root PATH` | `SCIENCE_ROOT` | Directory path | Root searched recursively for Science FITS/FZ archives. |
| `--dq-root PATH` | `DQ_ROOT` | Directory path | Root searched recursively for matching DQ FITS/FZ archives. |
| `--output-root PATH` | `OUTPUT_ROOT` | Directory path | Parent of target directories and generated exposure lists. |
| `--dataset TARGET:PREFIX` | `DATASETS` | Repeatable pair with exactly one `:` | Adds one dataset. `TARGET` must be one unique directory name, not `.`, `..`, or a path; `PREFIX` must be non-empty. Do not mix with `--target` or `--prefix`. |
| `--target NAME` | First configured dataset target | One directory name | Legacy single-dataset override. Do not mix with `--dataset`. |
| `--prefix TEXT` | First configured dataset prefix | Non-empty text | Legacy single-dataset archive prefix. Do not mix with `--dataset`. |
| `--contains TEXT` | `CONTAINS` | Repeatable non-empty text | Case-sensitive archive-basename substring filters with OR semantics. An empty configured list disables token filtering. |
| `--existing MODE` | `fail` | `fail`, `resume`, `overwrite` | Initializer publication policy: reject an existing target, reuse verified completed output, or replace output through staging. |
| `--f77-max-path N` | Standard `150` | Non-negative integer | Reject generated paths longer than the legacy reader can safely consume. `0` disables the check. |
| `--expo-list PATH` | `EXPO_LIST=""` | File path | Exposure list for main-only mode. With multiple main-only datasets, one explicit list is rejected; omit it to derive `OUTPUT_ROOT/expo_<target>.list` for each dataset. |
| `LEGACY_EXPO_LIST` | None | One positional path | Compatibility alias for `--expo-list`. Only one positional value is accepted, and explicit `--expo-list` wins. |
| `--help` | Off | No value | Print usage and exit successfully without running a phase. |

At least one dataset is required when init or main runs. Extcat-only mode does
not require a dataset. In chained init/main mode, each successful
`process_init` result overrides any external or configured exposure-list path
for that dataset.

## `ProcessConfig.hpp` workflow defaults

These values seed `RuntimeOptions`. Except for the noted Standard/Lite
difference, both variants expose the same fields.

### Phase and extcat defaults

| Parameter | Standard default | Lite default | Purpose / valid values |
|:---|:---|:---|:---|
| `RUN_PROCESS_EXTCAT` | `false` | `false` | Default `process_extcat` phase switch. |
| `RUN_PROCESS_INIT` | `false` | `false` | Default `process_init` phase switch. |
| `RUN_PROCESS_MAIN` | `true` | `true` | Default `process_main` phase switch. |
| `EXTCAT_INPUT_DIRECTORY` | Empty | Empty | Raw catalog root; must be set if extcat is enabled. |
| `EXTCAT_OUTPUT_DIRECTORY` | Reference to `LensingConfig::SOURCE_CAT` | Same | Read-only alias used to seed `RuntimeOptions::extcat_output_directory`; edit `SOURCE_CAT`, not this alias. |
| `EXTCAT_FILENAME_TOKENS` | Empty | Empty | Basename substring filters; empty accepts all files. |
| `EXTCAT_RECURSIVE` | `true` | `true` | Recursive discovery switch. |
| `EXTCAT_DELIMITER` | `"auto"` | Same | `auto`, `whitespace`, `comma`, or `tab`. |
| `EXTCAT_HEADER_MODE` | `"auto"` | Same | `auto`, `present`, or `absent`. |
| `EXTCAT_MALFORMED_POLICY` | `"fail"` | Same | `fail` or `skip`. |
| `EXTCAT_EXISTING_POLICY` | `"fail"` | Same | `fail` or `overwrite`. |
| `EXTCAT_CHUNK_MIB` | `64` | `64` | Positive MPI task size in MiB. |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `false` | `false` | Enables `EXTCAT_INPUT_COLUMNS_ONE_BASED`. When `false`, all raw fields retain their input order. |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `{1, ..., 18}` | Same | Any non-empty ordered vector of positive raw-field indices, used only when explicit selection is enabled. Its length sets output width and repeated values repeat fields. |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` | `false` | `false` | When `false`, unique `ra`/`dec` header names are preferred and the configured indices are fallback values. When `true`, always use the configured raw indices. |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `1` | `5` | Positive raw RA fallback used for sky tiling; must differ from Dec. A recognized `ra` header takes precedence unless explicit coordinate mode is enabled. |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `2` | `6` | Positive raw Dec fallback used for sky tiling; must differ from RA. A recognized `dec` header takes precedence unless explicit coordinate mode is enabled. |

### Initializer and exposure-list defaults

| Parameter | Standard default | Lite default | Purpose / valid values |
|:---|:---|:---|:---|
| `SCIENCE_ROOT` | `/lustre/home/acct-phyzj/share/DES/g` | Same | Default Science archive root. |
| `DQ_ROOT` | `/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask` | Same | Default DQ archive root. |
| `OUTPUT_ROOT` | `/lustre/home/acct-phyzj/share/DES/g_band_v1` | Same | Default output/list root. |
| `DATASETS` | `{{"gband", "c4d_"}}` | `{{"g2019", "c4d_19"}}` | Vector of `DatasetSpec{target, prefix}` pairs. Targets must be unique safe directory names and prefixes non-empty. |
| `CONTAINS` | `{"v1"}` | Same | OR-matched Science/DQ archive basename tokens. Empty means no token filter. |
| `EXISTING` | `"fail"` | Same | `fail`, `resume`, or `overwrite`. |
| `F77_MAX_PATH` | `150` | Non-negative generated-path limit; `0` disables it. |
| `EXPO_LIST` | Empty | Empty | Optional default main-only exposure-list path. |

### `RuntimeOptions` adapter fields

`DatasetSpec::target` and `DatasetSpec::prefix` are the two fields represented
by `DATASETS` and `--dataset`. `RuntimeOptions` mirrors the public defaults as
shown below; this is useful when calling `process_extcat`, `process_init`, or
`process_main` from another C++ driver.

| `RuntimeOptions` field | Seeded from / runtime override |
|:---|:---|
| `run_process_extcat` | `RUN_PROCESS_EXTCAT` / `--run-extcat` |
| `run_process_init` | `RUN_PROCESS_INIT` / `--run-init` |
| `run_process_main` | `RUN_PROCESS_MAIN` / `--run-main` |
| `extcat_input_directory` | `EXTCAT_INPUT_DIRECTORY` / `--extcat-input` |
| `extcat_output_directory` | `SOURCE_CAT` through `EXTCAT_OUTPUT_DIRECTORY` / `--extcat-output` |
| `extcat_filename_tokens` | `EXTCAT_FILENAME_TOKENS` / repeatable `--extcat-contains` |
| `extcat_recursive` | `EXTCAT_RECURSIVE` / `--extcat-recursive` |
| `extcat_delimiter` | `EXTCAT_DELIMITER` / `--extcat-delimiter` |
| `extcat_header_mode` | `EXTCAT_HEADER_MODE` / `--extcat-header` |
| `extcat_malformed_policy` | `EXTCAT_MALFORMED_POLICY` / `--extcat-malformed` |
| `extcat_existing_policy` | `EXTCAT_EXISTING_POLICY` / `--extcat-existing` |
| `extcat_chunk_mib` | `EXTCAT_CHUNK_MIB` / `--extcat-chunk-mib` |
| `extcat_use_explicit_columns` | `EXTCAT_USE_EXPLICIT_COLUMNS`; set to `true` by `--extcat-columns` |
| `extcat_input_columns_one_based` | `EXTCAT_INPUT_COLUMNS_ONE_BASED` / `--extcat-columns` |
| `extcat_use_explicit_coordinate_columns` | `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS`; set to `true` by either coordinate-column CLI option |
| `extcat_ra_column_one_based` | `EXTCAT_RA_COLUMN_ONE_BASED` / `--extcat-ra-column` |
| `extcat_dec_column_one_based` | `EXTCAT_DEC_COLUMN_ONE_BASED` / `--extcat-dec-column` |
| `science_root` | `SCIENCE_ROOT` / `--science-root` |
| `dq_root` | `DQ_ROOT` / `--dq-root` |
| `output_root` | `OUTPUT_ROOT` / `--output-root` |
| `datasets` | `DATASETS` / `--dataset` or legacy `--target` and `--prefix` |
| `contains` | `CONTAINS` / repeatable `--contains` |
| `existing` | `EXISTING` / `--existing` |
| `f77_max_path` | `F77_MAX_PATH` / `--f77-max-path` |
| `expo_list` | `EXPO_LIST` / `--expo-list` or positional `LEGACY_EXPO_LIST` |
| `external_expo_list_supplied` | Parser state recording an explicit or positional exposure list; not an independent user parameter. |
| `help_requested` | Parser state set by `--help`; not an independent user parameter. |

## `LensingConfig.hpp` scientific parameters

All parameters in this section are compile-time values. Rebuild after editing.
“Derived” values should normally follow their source parameter instead of being
tuned independently. “Reserved” values are retained for source/F77
compatibility but have no active effect in the current C++ implementation.

### Mathematical constants, dimensions, and paths

| Parameter | Standard / Lite default | Status | Function and constraints |
|:---|:---|:---|:---|
| `pi` | `3.14159265358979323846` | Derived mathematical constant | Value of π. |
| `arc_convert` | `pi / 180` | Derived | Degrees-to-radians conversion. |
| `npx` | `3000` | Active | Maximum/fallback image width and allocation width for very-local PSF maps; must exceed real chip width. |
| `npy` | `5000` | Active | Maximum/fallback image height and allocation height for very-local PSF maps; must exceed real chip height. |
| `strl` | `150` | Reserved | Legacy fixed string-length compatibility constant; current C++ path handling uses `std::string`. |
| `ASTROMETRY_CAT` | `/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted` | Active path | Gaia tile root used when nontrivial astrometry is selected. |
| `SOURCE_CAT` | `/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat` | Primary path with runtime override | One-degree external catalog directory. It seeds `EXTCAT_OUTPUT_DIRECTORY`; `--extcat-output` writes the effective value back before processing starts. |
| `FLAT_PATH` | `/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014` | Standard only | Root for legacy flat/mask FITS files when their branches are selected. |
| `PSF_PATH` | `"hahahaha"` | Standard only; placeholder default | Directory containing `PSF.fits` when `ext_PSF=1`; must be replaced before using that branch. |

### Stage and branch selection

`PROCESS_stage` is the product of the prime numbers assigned to enabled stages.
A stage runs when `PROCESS_stage % prime == 0`.

| Prime | Stage | C++ module |
|---:|:---|:---|
| 2 | Preprocessing | `PreProcess` |
| 3 | Astrometric calibration | `Astrometry` |
| 5 | Source detection/extraction | `SourceExtractor` |
| 7 | First Fourier transform | `FourierTransformSt1` |
| 11 | PSF modeling | `PSFModel` and optional `PSFRecons` |
| 13 | Second Fourier transform | `FourierTransformSt2` |
| 17 | Shear measurement | `ShearMeasurement` |
| 19 | Exposure information | `ExposureInfo` |
| 23 | Catalog combination | `CatalogCombiner` |

| Parameter | Standard default / options | Lite behavior | Function and dependencies |
|:---|:---|:---|:---|
| `PROCESS_stage` | Product of all nine primes | Same and editable | Enabled-stage product. Omitting a prime disables that stage, but required intermediates must already exist if dependent later stages remain enabled. |
| `ASTROMETRY_trivial` | `0`; `0` Gaia fit, `1` trivial WCS | Frozen to `0`, constant removed | Selects astrometric calibration branch. `0` requires `ASTROMETRY_CAT`. |
| `include_FLAT` | `0`; `0` off, `1` apply flat | Frozen to `0`, constant removed | Enables legacy super-flat multiplication. Requires valid `FLAT_PATH` products. |
| `include_Mask` | `2`; `0` none, `1` legacy mask, `2` per-chip DQ, `3` both | Frozen to `2`, constant removed | Selects mask sources. Values `2/3` require initializer-compatible DQ masks. |
| `ext_cat` | `1`; `0` internal detection, `1` external catalog | Frozen to `1`, constant removed | Selects source positions. Value `1` requires `SOURCE_CAT` and the expected post-RA schema. |
| `ext_PSF` | `0`; `0` derive from frame stars, `1` external fixed PSF | Frozen to `0`, constant removed | Value `1` skips star FFT/modeling and requires `PSF_PATH/PSF.fits`. |
| `deblending` | `1`; `0` off, `1` on | Frozen to `1`, constant removed | Applies redshift-consistency deblending to external-catalog sources. |
| `PSF_type` | `1`; `1` local polynomial, `2` very-local map/interpolation | Frozen to `1`, constant removed | Chooses PSF model representation. Type `2` uses `step_psf` and `n_neighbor`. |
| `PSF_Ms` | `0`; `0` local model only, `1` add multi-scale/PCA residual reconstruction | Frozen to `0`, constant removed | PCA mode is available only in Standard and requires its Stage 5 products and PCA constants. |

### External catalogs, background, and PSF fitting

| Parameter | Standard default | Lite default | Function and options |
|:---|---:|---:|:---|
| `ext_cat_columns_before_ra` | `4` | `4` | Number of whitespace-delimited fields skipped before RA in each `SOURCE_CAT` row. Use `4` for the DES flag columns or `0` when RA is first. The reader then expects Dec, five magnitude/error pairs, `dnf_z`, and `dnf_zsigma`; arbitrary extcat projection is compatible with `process_main` only when it preserves this layout. |
| `CCD_split` | `2` | `2` | Number of horizontal amplifier regions used for background/noise fits. The current storage supports `1` (whole chip) or `2` (two amplifiers); other values are unsafe. |
| `blocksize` | `200` | `200` | Pixel width/height of blocks sampled for background estimation. |
| `nct` | `12` | `12` | Number of rectangular-monomial terms in the background surface fit. Must not exceed the available stable background samples. |
| `ncx` | `3` | `3` | Number of successive x powers in the background basis before y power increments. With `nct=12`, the basis spans `x^0..2` for `y^0..3`. |
| `psf_order` | `8` | `8` | Reserved legacy PSF-order constant; currently unused. |
| `npo` | `64` | `64` | Reserved legacy PSF selection size; currently used only to derive `nstar_min`. |
| `npox` | `8` | `8` | Reserved legacy PSF basis-width constant; currently unused. |
| `nstar_min` | `npo * 3 / 2 = 96` | Same | Base exposure-wide star threshold. PSF star selection rejects an exposure when the total candidate count is below `2 * nstar_min` (192 by default). |
| `npl` | `10` | `10` | Number of ordered 2D polynomial terms fitted per PSF Fourier pixel; `10` includes terms through total degree 3. |
| `nplx` | `2` | `2` | Compatibility value copied into PSF routines but currently not consumed downstream. |
| `nstar_min_local` | `16` | `16` | Minimum finite stars required for one local chip PSF fit. |
| `step_psf` | `100` | Not present | Standard very-local PSF map grid spacing in pixels; used only with `PSF_type=2`. |
| `n_neighbor` | `5` | Not present | Standard number of nearest stars used by the very-local PSF branch. |

### Stamp geometry, detection, and Fourier power

| Parameter | Standard default | Lite default | Function and constraints |
|:---|---:|---:|:---|
| `ns` | `64` | `64` | Square source/star stamp width in pixels; FFT routines assume the configured dimensions consistently. |
| `nsns` | `ns * ns = 4096` | Same | Derived stamp pixel count. |
| `chip_margin` | `8` | `8` | Extra pixels around the half-stamp extraction region. |
| `ns_2` | `ns / 2 = 32` | Same | Derived half-stamp size. |
| `nl_2` | `ns_2 + chip_margin = 40` | Same | Derived half-width of the larger extraction region. |
| `nl` | `2 * nl_2 = 80` | Same | Derived larger extraction width. |
| `flag_thresh` | `3` | `3` | Maximum reserved/flagged edge distance used when accepting a source stamp. |
| `chip_edge_margin` | `chip_margin = 8` | Same | Derived alias used to reject sources too near chip edges. |
| `dz_thresh` | `0.1` | `0.1` | Maximum redshift difference for keeping neighboring entries during external-catalog deblending. |
| `source_thresh` | `2.0` | `2.0` | Detection/defect connected-pixel threshold in noise-sigma units. |
| `core_thresh` | `4.0` | `4.0` | Detection-core/peak threshold in noise-sigma units. |
| `flat_thresh` | `0.01` | Not present | Reserved Standard legacy flat threshold; currently unused. |
| `area_max` | `ns * ns = 4096` | Same | Derived maximum connected-region workspace/area. |
| `area_thresh` | `6` | `6` | Minimum connected-region pixel count. |
| `gal_smooth` | `0` | `2` | Galaxy/noise FFT power smoothing: `0` none, `1` linear 5x5 hole-aware, `2` log-power 5x5 hole-aware. This intentional default difference changes main-galaxy processing. |
| `star_smooth` | `2` | `2` | Star FFT power smoothing with the same `0/1/2` modes. In star power normalization, `0` uses four neighboring central pixels while values `>=1` use the central pixel. |
| `SNR_PSF` | `100.0` | `100.0` | S/N threshold used to select PSF-star candidates; some preliminary selection uses half this value. |
| `saturation_thresh` | `25000.0` | `25000.0` | Raw-pixel saturation cutoff and normalized peak rejection reference. |
| `pixel_size` | `0.2628` arcsec | Same | Detector pixel scale used to convert PSF sizes to angular units. |

### Catalog and memory dimensions

| Parameter | Default | Status | Function and constraints |
|:---|---:|:---|:---|
| `len_g` | `40` | Active | Galaxy stamps per FITS layout row/block used by stamp I/O. |
| `len_s` | `15` | Active | Star stamps per FITS layout row/block. |
| `ngal_max` | `4000` | Active capacity | Maximum galaxies stored per chip. Larger detected catalogs are truncated at this limit. |
| `nstar_max` | `2000` | Active capacity | Maximum stars stored per chip. Memory use includes arrays scaling as `nstar_max^2`. |
| `npara` | `25` | Active layout | Number of per-source/per-star parameter slots in internal tables. Must cover all configured indices. |
| `len_sam` | `50` | Active | Number of exposure-wide selected-star stamps placed in one FITS layout row/block. |
| `npd` | `33` | Active layout | Number of astrometric PU distortion coefficients per coordinate. Must match astrometry file serialization and fitting code. |
| `NMAX_EXPO` | `25000` | Active capacity | Maximum exposure records allocated for aggregation. |
| `NMAX_CHIP` | `62` | Active capacity | Maximum chips represented in fixed-capacity PSF/exposure arrays. |

### Mode-bar noise-plane estimator

These F6 parameters act together. They should normally remain synchronized with
the validated estimator convention rather than be tuned independently.

| Parameter | Default | Function |
|:---|---:|:---|
| `sig_blocksize` | `200` | Pixel side length of one robust seed block. |
| `sig_block_max` | `sig_blocksize * sig_blocksize = 40000` | Derived maximum pixels stored per block. |
| `sig_max_blocks` | `2048` | Maximum block seeds retained per amplifier. |
| `sig_min_block_pixels` | `1000` | Minimum valid pixels required for one block seed. |
| `sig_min_block_triples` | `1000` | Minimum valid pixel triples required for a block estimate. |
| `sig_min_blocks` | `4` | Minimum accepted blocks required for the estimator. |
| `sig_hist_nbin` | `256` | Histogram bin count for locating the sky mode. |
| `sig_hist_range` | `6.0` | Histogram half/range scale in robust-width units. |
| `sig_min_mode_count` | `500` | Minimum sample count supporting mode estimation. |
| `sig_min_lower_count` | `1000` | Minimum lower-side sample count for width estimation. |
| `sig_lower_quantile` | `0.3173105` | Lower-side quantile used by the mode-width estimator. |
| `sig_clip_k` | `3.0` | Symmetric clipping threshold in sigma units. |
| `sig_rdil` | `2` | Radius used to dilate the estimator's private brightness mask. |
| `sig_clip_niter` | `2` | Number of clipped plane-fit iterations. |
| `sig_min_fit_triples` | `1000` | Minimum triples entering a noise-plane fit. |
| `sig_min_fit_frac` | `0.20` | Minimum surviving fraction of possible fit triples. |
| `sig_median_ratio` | `1.2678405` | Calibration from robust lower-side width to sigma convention. |
| `sig_plane_min` | `1.0e-8` | Minimum positive fitted noise-plane value. |
| `sig_max_plane_ratio` | `4.0` | Maximum allowed variation ratio across a fitted plane. |
| `sig_pivot_min` | `1.0e-8` | Minimum numerical pivot accepted by the small plane solve. |
| `sig_scale_s1` | `0.673475` | Retained Stage-1 calibration candidate. |
| `sig_scale_s2` | `1.027786` | Stage-2 calibration candidate used by the current pipeline. |
| `sig_scale` | `sig_scale_s2` | Active derived selector converting the fitted plane to the published `2*sigma^2` convention. |

### Internal catalog column indices

These are zero-based positions in the C++ per-source result rows. They are not
the 18 raw external-catalog projection indices. Changing them changes the
internal/output layout and requires coordinated reader/writer changes.

| Parameter | Index | Field / status | Parameter | Index | Field / status |
|:---|---:|:---|:---|---:|:---|
| `isig` | 3 | Reserved; currently unused | `istar` | 4 | Star indicator / classification field |
| `ipeak` | 4 | Reserved alias; currently unused | `i_imax` | 5 | Peak x index |
| `i_jmax` | 6 | Peak y index | `ih_flux` | 7 | Half-light/source flux field |
| `ih_area` | 8 | Source area field | `iflag` | 9 | Quality flag |
| `iPSF` | 10 | PSF/SNR-related stored field | `iSNR_F` | 11 | Fourier S/N field |
| `ira` | 12 | Right ascension | `idec` | 13 | Declination |
| `igf1` | 14 | Field-distortion component 1 | `igf2` | 15 | Field-distortion component 2 |
| `ig1` | 16 | Fourier_Quad shear estimator 1 | `ig2` | 17 | Fourier_Quad shear estimator 2 |
| `ide` | 18 | Fourier_Quad normalization estimator | `ih1` | 19 | Higher-order estimator 1 |
| `ih2` | 20 | Higher-order estimator 2 | `icos2` | 21 | Spin-2 cosine term |
| `isin2` | 22 | Spin-2 sine term | `iparity` | 23 | WCS parity |

### Calibration and camera geometry

| Parameter | Default | Function and constraints |
|:---|---:|:---|
| `g1_c` | `-0.001` | Additive correction applied to field-distortion component 1 during final catalog combination. |
| `g2_c` | `-0.0003` | Additive correction applied to field-distortion component 2 during final catalog combination. |
| `chi2_thresh` | `0.01` | Maximum accepted exposure/PSF chi-square diagnostic during catalog combination. |
| `chipnx` | `2046` | Nominal science CCD width used to normalize PSF coordinates. Must match the instrument geometry. |
| `chipny` | `4094` | Nominal science CCD height used to normalize PSF coordinates. Must match the instrument geometry. |
| `Camera_ccd_num` | `62` | Instrument CCD count used by exposure and PCA allocations. |

### Standard-only multi-scale/PCA PSF parameters

These values are compiled only by `cpp_Standard` and are active only when
`PSF_Ms=1`. `cpp_Lite` removes the PCA implementation and all of these
parameters.

| Parameter | Default | Function and constraints |
|:---|---:|:---|
| `rescale_size` | `1.2` | Target PSF size used to derive the residual rescaling factor. |
| `procs_pn` | `40` | Process-group size passed to PCA reconstruction scheduling. |
| `work_pn` | `10` | Concurrent worker count within the PCA scheduler grouping. |
| `nblocks` | `2` | Number of spatial blocks per CCD axis; default creates a 2x2 block grid. |
| `n_pcs` | `100` | Maximum number of residual principal components stored and fitted. |
| `npp6th` | `28` | Number of ordered 2D sixth-degree polynomial terms used for PCA coefficient surfaces. |
| `pca_negative_eigenvalue_threshold` | `-1.0e-5` | Eigenvalue below which a PCA covariance result is classified as invalid. |
| `nmax_star_pchip` | `1000000` | Reserved legacy per-chip PCA star capacity; currently unused. |

## Standard and Lite summary

The runtime command line is intentionally identical in both variants. Their
current configuration differences are:

- Standard default dataset: `gband:c4d_`; Lite: `g2019:c4d_19`.
- Standard `F77_MAX_PATH=150`; Lite `149`.
- Standard `gal_smooth=0`; Lite `gal_smooth=2`.
- Lite freezes Gaia astrometry, no flat, DQ-mask mode, external source catalogs,
  frame-derived PSF stars, deblending, local-polynomial PSF, and no PCA. The
  unselected Standard branches and their exclusive constants are absent.
- Standard retains the optional external PSF, very-local PSF, and PCA branches.

When moving a parameter file between variants, preserve these intentional
differences instead of copying one `LensingConfig.hpp` over the other.
