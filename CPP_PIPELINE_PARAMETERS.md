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

### Stage-5 PSF star selection

These controls are compile-time settings shared by Standard and Lite:

| Symbol | Default | Meaning |
|---|---:|---|
| `psf_exposure_min_candidates` | `60` | Minimum quality-valid candidates before exposure selection. |
| `psf_fwhm_hist_bins` | `128` | FWHM-locus histogram bins. |
| `psf_fwhm_locus_sigma` | `4.0` | Robust FWHM-locus width multiplier. |
| `psf_fwhm_locus_min_samples` | `30` | Minimum samples for a valid exposure locus. |
| `psf_minchi_reference_fraction` | `1/3` | Largest exposure-wide locus fraction eligible as threshold references. |
| `psf_minchi_reference_max_per_chip` | `5` | Maximum references retained per chip within that top fraction. |
| `psf_minchi_sigma_cut` | `4.0` | Upper-tail cut for the reference-all pair sample. |
| `PsfGroupingType` | `2` | `1` legacy threshold graph; `2` survivor-only mutual KNN. |
| `psf_knn_k` | `8` | Exact neighbours retained in mutual-KNN mode. |
| `psf_group_merge_ratio` | `0.30` | Secondary/main group size ratio. |
| `psf_group_merge_min_gaia` | `2` | Gaia matches required for a secondary group. |
| `psf_press_rejection_enabled` | `true` | Enable optional standardized-PRESS removal/refit. |
| `psf_press_sigma_cut` | `4.0` | Exposure cut for leverage-standardized PRESS. |
| `psf_press_max_removals` | `5` | Maximum proposed removals allowed per chip. |
| `psf_loo_min_denom` | `1e-6` | Minimum accepted analytic-LOO denominator `1-h`. |

Every same-chip locus pair contributes to the endpoint `min_chi` values. Only
unique pairs with at least one capped large-size reference contribute to the
common minChi threshold. The Type-1 graph continues to estimate its private
legacy all-FWHM-pair threshold.

Raw PRESS is the central-window analytic-LOO RMS. Rejection uses
`raw_press * sqrt(1 - leverage)` without another amplitude normalization. The
first fit remains final when rejection is disabled, no outlier is found, more
than five are flagged, fewer than 16 would remain, or the temporary refit
fails. Thus `psf_press_rejection_enabled=false` disables only rejection, not
the first fit, leverage, analytic LOO, raw/standardized scores, or final LOO
residual generation.

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
