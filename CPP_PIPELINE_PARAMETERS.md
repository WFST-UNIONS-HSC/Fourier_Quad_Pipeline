# C++ Pipeline Parameter Reference

This file describes the C++ programs in this repository. They use compiled
defaults plus CLI overrides; there is no INI configuration layer.

## Run-time CLI

Use `./Fourier_Quad_Pipe --help` for the authoritative current list. Options
accept `--name value` and `--name=value`; booleans accept `true/false`, `1/0`,
and `on/off`.

### Workflow and downstream paths

| CLI | Standard default | Lite default | Meaning |
|---|---:|---:|---|
| `--run-extcat` | `false` | `false` | External-catalog tiling. |
| `--run-init` | `true` | `true` | Archive initialization. |
| `--run-main` | `true` | `true` | Nine-stage numerical pipeline. |
| `--run-rearr` | `true` | `false` | Spatial catalog rearrangement. |
| `--run-fd` | `true` | `false` | FD shear test. |
| `--expo-list` | empty | empty | Downstream-only exposure list. |
| `--rearr-output-dir` | `baked` | `baked` | Rearrangement directory name. |
| `--rearr-output-base` | empty | empty | Empty means dataset root. |
| `--rearr-list-name` | `cat_gband_ori.list` | same | Rearranged list filename. |
| `--rearr-list-dir` | empty | empty | Empty means input-list parent. |
| `--fd-expo-list` | empty | empty | Optional FD list override. |
| `--fd-output-dir` | `fdout` | `fdout` | FD result directory name. |
| `--fd-output-base` | empty | empty | Empty means dataset root. |

At least one phase must be enabled.

### External catalog

| CLI | Default | Meaning |
|---|---|---|
| `--extcat-input` | empty | Raw catalog root. |
| `--extcat-output` | compiled `SOURCE_CAT` | Tile output and effective pipeline catalog path. |
| `--extcat-contains` | no filter | Repeatable basename substring; OR semantics. |
| `--extcat-recursive` | `true` | Scan subdirectories. |
| `--extcat-delimiter` | `auto` | `auto`, `whitespace`, `comma`, or `tab`. |
| `--extcat-header` | `auto` | `auto`, `present`, or `absent`. |
| `--extcat-columns` | pass-through | Ordered one-based projection. |
| `--extcat-ra-column` | `5` | Raw one-based RA position. |
| `--extcat-dec-column` | `6` | Raw one-based Dec position. |
| `--extcat-zp-column` | `17` | Raw one-based photo-z position. |
| `--extcat-chunk-mib` | `64` | Nominal MPI task size. |
| `--extcat-malformed` | `fail` | `fail` or `skip`. |
| `--extcat-existing` | `fail` | `fail` or `overwrite`. |

With explicit projection, RA, Dec, and photo-z must be present. The output
directory cannot equal or be nested below the input directory.

### Initializer

| CLI | Default source | Meaning |
|---|---|---|
| `--science-root` | `config/InitConfig.hpp` | Science archive root. |
| `--dq-root` | `config/InitConfig.hpp` | DQ archive root. |
| `--output-root` | `config/InitConfig.hpp` | Writable processing root. |
| `--dataset TARGET:PREFIX` | compiled list | Repeatable dataset; runs are sequential. |
| `--target`, `--prefix` | none | Legacy single-dataset form; do not mix with `--dataset`. |
| `--contains` | compiled list | Repeatable archive basename token; OR semantics. |
| `--existing` | `fail` | `fail`, `resume`, or `overwrite`. |
| `--f77-max-path` | `150` | Generated-path limit; zero disables it. |

The first explicit dataset/filter occurrence replaces the corresponding
compiled list; later occurrences append.

## Compile-time settings

CLI only overrides fields represented in `ProcessConfig::RuntimeOptions`.
Changing the settings below requires rebuilding the selected variant.

| Header | Important controls |
|---|---|
| `config/LensingConfig.hpp` | `PROCESS_stage`, Standard branch switches, astrometry/source/flat/PSF paths, stamp geometry, detection thresholds, PSF grouping, noise construction, catalog indices, calibration, CCD geometry, and Standard PCA settings |
| `config/ProcessRearrConfig.hpp` | 0.1-degree grid, target partition size, output names/precision, missing/malformed policies |
| `config/FDConfig.hpp` | FD/PDF/jackknife mode, bins, catalog cuts, stellar-locus selection, excluded CCDs, and chip-edge mask |

Common values include:

| Symbol | Default | Meaning |
|---|---:|---|
| `PROCESS_stage` | `223092870` | All nine prime-gated stages. |
| `ns` | `64` | Stamp/Fourier side length. |
| `NstampType` | `1` | Blank-noise stamp (`1`) or local covariance power (`2`). |
| `source_thresh`, `core_thresh` | `2.0`, `4.0` | Detection thresholds. |
| `SNR_PSF` | `100` | PSF-star S/N threshold. |
| `ngal_max`, `nstar_max` | `4000`, `2000` | Initial vector reservation hints, not truncation limits. |
| `pixel_size` | `0.2628` | Arcsec per pixel. |
| `SKY_GRID_DEGREES` | `0.1` | Rearrangement grid size. |
| `TARGET_SUBCAT_ROWS` | `500000` | Target rows per spatial partition. |
| `fd_num`, `PDF_BINS` | `21`, `4` | FD outer and inner bin counts. |

Standard retains eight selectable branches. Lite fixes them to Gaia
astrometry, no flat, DQ mask mode 2, external catalog, frame-star PSF,
deblending, local PSF, and no PCA. Missing Lite branches cannot be restored by
adding their constants.

## Output catalog

Stage 7 writes 24 fields in this order:

| # | Field | # | Field |
|---:|---|---:|---|
| 1 | PSF polynomial chi-square | 13 | RA |
| 2 | source x | 14 | Dec |
| 3 | source y | 15 | field distortion `g1` |
| 4 | local noise sigma | 16 | field distortion `g2` |
| 5 | available PSF-star count | 17 | Fourier_Quad `g1` |
| 6 | peak x | 18 | Fourier_Quad `g2` |
| 7 | peak y | 19 | response `de` |
| 8 | half-light flux | 20 | higher-order `h1` |
| 9 | source area | 21 | higher-order `h2` |
| 10 | quality flag | 22 | spin-2 cosine |
| 11 | local PSF size | 23 | spin-2 sine |
| 12 | Fourier S/N | 24 | WCS parity |

Stage 9 appends exposure chi-square as field 25. With the default external
schema, the final row is 18 external fields + one CCD number + 25 pipeline
fields = 44 fields. With projection, use `external_width + 1 + 25`.
