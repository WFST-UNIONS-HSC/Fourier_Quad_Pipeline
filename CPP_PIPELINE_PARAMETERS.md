# C++ Pipeline Parameter Reference

This reference follows the live C++ source. Each configuration namespace has one
complete table; Standard and Lite values are compared in the same row. The C++
driver seeds `ProcessConfig::RuntimeOptions` from these compiled defaults and
only the options named in the **CLI override** column can change a value without
rebuilding.

## Centralized path configuration

Each variant has its own `config/pathconfig.hpp`. It is the sole physical source
for fixed input/output paths, workflow list/output names, rearrangement filenames,
and fixed relative output-directory layouts. The established namespaces remain
unchanged, so existing call sites still use names such as
`LensingConfig::ASTROMETRY_CAT` and
`ProcessConfig::REARR_OUTPUT_DIRECTORY`.

| Namespace | Definitions physically owned by `pathconfig.hpp` | Coupling and runtime rule |
|---|---|---|
| `LensingConfig` | `ASTROMETRY_CAT`, `SOURCE_CAT_DEFAULT`; Standard only: `FLAT_PATH`, `PSF_PATH` | `ASTROMETRY_CAT`, `FLAT_PATH`, and `PSF_PATH` are compile-time only. `SOURCE_CAT_DEFAULT` seeds the runtime external-catalog directory, which `--extcat-output` may replace. |
| `AstroCatConfig` | `ASTROCAT_INPUT_DIRECTORY`, `ASTROCAT_OUTPUT_DIRECTORY` | The compiled output default is deliberately initialized as `LensingConfig::ASTROMETRY_CAT`. `--astrocat-output` changes only the runtime producer destination and does not update the Stage-1 consumer path. |
| `ExtCatConfig` | `EXTCAT_INPUT_DIRECTORY`, `EXTCAT_OUTPUT_DIRECTORY` | The compiled output default remains `LensingConfig::SOURCE_CAT_DEFAULT`. The runtime catalog directory is shared by `process_extcat` output and `process_main` input. |
| `InitConfig` | `SCIENCE_ROOT`, `DQ_ROOT`, `OUTPUT_ROOT` | These seed `RuntimeOptions` and have CLI overrides. |
| `ProcessConfig` | `EXPO_LIST`, `REARR_OUTPUT_DIRECTORY`, `REARR_OUTPUT_BASE_DIRECTORY`, `REARRANGED_EXPO_LIST_FILENAME`, `REARRANGED_EXPO_LIST_DIRECTORY`, `FD_EXPO_LIST`, `FD_OUTPUT_DIRECTORY`, `FD_OUTPUT_BASE_DIRECTORY` | These seed `RuntimeOptions` and have CLI overrides. |
| `ProcessRearrConfig` | `SKIP_DIRECTORY_NAME`, `SUBCAT_PREFIX`, `SUBCAT_EXTENSION`, `SUMMARY_FILENAME` | No runtime override; edit the selected variant and rebuild. |
| `OutputLayout` | `NON_CHIP_BASE_DIRECTORIES`, `CHIP_PRODUCT_DIRECTORIES` | No runtime override; these are fixed relative directory contracts used by initialization and processing. |

Editing `pathconfig.hpp` changes compiled defaults and therefore requires a clean
rebuild. A CLI override changes only the corresponding `RuntimeOptions` copy; it
does not mutate the header constants or rewrite their compiled relationships.

`N/A — removed in Lite` means the Lite source physically removed the alternate
branch. Adding the constant back does not restore that behavior. A **derived
parameter** is retained for completeness but must not be edited directly; change
its source parameter and keep the associated assertions and consumers consistent.

## `ProcessConfig` (`config/ProcessConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `RUN_PROCESS_ASTROCAT` | `bool` | `false` | `false` | `--run-astrocat` | Boolean | Run Gaia-catalog tiling. | Select phases per run with CLI. | CLI override; rebuild not required |
| `RUN_PROCESS_EXTCAT` | `bool` | `false` | `false` | `--run-extcat` | Boolean | Run external-catalog tiling. | Select phases per run with CLI. | CLI override; rebuild not required |
| `RUN_PROCESS_INIT` | `bool` | `true` | `true` | `--run-init` | Boolean | Run archive initialization. | Select phases per run with CLI. | CLI override; rebuild not required |
| `RUN_PROCESS_MAIN` | `bool` | `true` | `true` | `--run-main` | Boolean | Run the nine-stage numerical pipeline. | Select phases per run with CLI. | CLI override; rebuild not required |
| `RUN_PROCESS_REARR` | `bool` | `true` | `false` | `--run-rearr` | Boolean | Run spatial catalog rearrangement. | Select phases per run with CLI. | CLI override; rebuild not required |
| `RUN_PROCESS_FD` | `bool` | `true` | `false` | `--run-fd` | Boolean | Run the field-distortion shear test. | Select phases per run with CLI. | CLI override; rebuild not required |
| `EXPO_LIST` | `const char*` | empty | empty | `--expo-list` or one positional path | Exposure-list path | Default downstream exposure list. | Prefer CLI for each run. | CLI override; rebuild not required |
| `REARR_OUTPUT_DIRECTORY` | `const char*` | `"baked"` | `"baked"` | `--rearr-output-dir` | Directory name/path | Rearranged-catalog output directory. | Change when output layout changes. | CLI override; rebuild not required |
| `REARR_OUTPUT_BASE_DIRECTORY` | `const char*` | empty | empty | `--rearr-output-base` | Empty = dataset root | Optional rearrangement output base. | Change when outputs live outside the dataset root. | CLI override; rebuild not required |
| `REARRANGED_EXPO_LIST_FILENAME` | `const char*` | `"cat_gband_ori.list"` | same | `--rearr-list-name` | Filename | Published rearranged exposure-list name. | Change for another naming convention. | CLI override; rebuild not required |
| `REARRANGED_EXPO_LIST_DIRECTORY` | `const char*` | empty | empty | `--rearr-list-dir` | Empty = input-list parent | Published rearranged exposure-list directory. | Change for another list location. | CLI override; rebuild not required |
| `FD_EXPO_LIST` | `const char*` | empty | empty | `--fd-expo-list` | Exposure-list path | Optional FD-specific input list. | Override when FD must use a different list. | CLI override; rebuild not required |
| `FD_OUTPUT_DIRECTORY` | `const char*` | `"fdout"` | same | `--fd-output-dir` | Directory name/path | FD result directory. | Change when output layout changes. | CLI override; rebuild not required |
| `FD_OUTPUT_BASE_DIRECTORY` | `const char*` | empty | empty | `--fd-output-base` | Empty = dataset root | Optional FD output base. | Change when outputs live outside the dataset root. | CLI override; rebuild not required |

The `WorkflowOptions`, `PipelineOptions`, `CatalogOptions`, `AstroCatOptions`,
`ExtCatOptions`, `InitOptions`, `RearrOptions`, `FDOptions`, and `RuntimeOptions` structs are the
runtime copies of values listed in this and the next two header tables. Their
members are parser state, not a second set of user defaults. `help_requested`
and `external_exposure_list_supplied` are internal parser flags.

## `InitConfig` (`config/InitConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `DatasetSpec::target` | `std::string` | `"gband"` in `DATASETS` | same | `--dataset TARGET:PREFIX` or `--target` | One non-empty directory name | Dataset output directory name. | Change for each dataset. | CLI override; rebuild not required |
| `DatasetSpec::prefix` | `std::string` | `"c4d_"` in `DATASETS` | same | `--dataset TARGET:PREFIX` or `--prefix` | Non-empty archive basename prefix | Select matching Science images. | Change for each dataset. | CLI override; rebuild not required |
| `SCIENCE_ROOT` | `const char*` | `/lustre/home/acct-phyzj/share/DES/g` | same | `--science-root` | Readable directory path | Science-image archive root. | Change for another site or archive. | CLI override; rebuild not required |
| `DQ_ROOT` | `const char*` | `/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask` | same | `--dq-root` | Readable directory path | DQ-mask archive root when the chosen configuration uses DQ masks. | Change for another site; irrelevant only when DQ access is disabled. | CLI override; rebuild not required |
| `OUTPUT_ROOT` | `const char*` | `/lustre/home/acct-phyzj/share/DES/g_band_v1` | same | `--output-root` | Writable directory path | Parent of dataset trees and generated exposure lists. | Change for every deployment. | CLI override; rebuild not required |
| `DATASETS` | `std::vector<DatasetSpec>` | `{{"gband", "c4d_"}}` | same | Repeatable `--dataset`; legacy `--target` + `--prefix` | Unique targets with non-empty prefixes | Datasets processed sequentially. | Change for another dataset set. | CLI override; rebuild not required |
| `CONTAINS` | `std::vector<std::string>` | `{"v1"}` | same | Repeatable `--contains` | Non-empty case-sensitive basename tokens; OR matching | Additional Science-image and DQ-mask archive filename filters. | Change when archive naming changes. | CLI override; rebuild not required |
| `EXISTING` | `const char*` | `"fail"` | same | `--existing` | `fail`, `resume`, `overwrite` | Existing-output policy for initialization. | Select intentionally per run. | CLI override; rebuild not required |
| `F77_MAX_PATH` | `int` | `150` | `150` | `--f77-max-path` | Non-negative; `0` disables | Compatibility guard for generated path lengths. | Change only for path-policy compatibility. | CLI override; rebuild not required |

## `AstroCatConfig` (`config/AstroCatConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `ASTROCAT_INPUT_DIRECTORY` | `const char*` | empty | empty | `--astrocat-input` | Readable flat directory | Raw Gaia files; each data row begins with RA and Dec. | Set when running `process_astrocat`. | CLI override; rebuild not required |
| `ASTROCAT_OUTPUT_DIRECTORY` | `std::string` | `LensingConfig::ASTROMETRY_CAT` | same | `--astrocat-output` | Writable directory that does not equal, contain, or sit below the input directory | Compiled default destination for one-degree Type-2 Gaia tiles. | Override independently for each publication. | No via CLI; yes after editing the header |
| `ASTROCAT_ADD_HEADER` | `bool` | `true` | `true` | `--astrocat-add-header` | `true` starts at the first line; `false` skips exactly one line per input file | Controls raw-input header handling; output tiles always contain `RA    DEC`. | Change to match the raw files. | CLI override; rebuild not required |
| `ASTROCAT_EXISTING_POLICY` | `const char*` | `"fail"` | same | `--astrocat-existing` | `fail`, `overwrite` | Existing generated-tile policy. | Select intentionally for reruns. | CLI override; rebuild not required |

In `pathconfig.hpp`, `ASTROCAT_OUTPUT_DIRECTORY` intentionally derives from
`LensingConfig::ASTROMETRY_CAT`; keep both symbols and that expression rather
than merging them. After `RuntimeOptions` is constructed, `--astrocat-output`
changes only the `process_astrocat` destination. It is not compared with or
propagated back to the compile-time Stage-1 consumer path.

The phase discovers only direct regular children of the input directory; it
does not recurse. It reads each complete file through dynamic MPI scheduling,
optionally skips exactly one first line, replaces commas with spaces, consumes
the first two parseable doubles, and silently skips rows without two doubles.
The input contract is finite sky coordinates with `0 <= RA <= 360` and
`-90 <= Dec <= 90`; exactly `RA=360` is stored as zero and exactly `Dec=90`
belongs to the last Dec tile. Exact and one-ULP duplicates in both coordinates
are removed, including duplicates across tile boundaries. Output files use
`des_y6_RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`, always begin with `RA    DEC`,
and contain round-trip-precision doubles. `overwrite` removes only files that
match this generated basename contract and preserves unrelated directory
content.

## `ExtCatConfig` (`config/ExtCatConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `EXTCAT_INPUT_DIRECTORY` | `const char*` | empty | empty | `--extcat-input` | Readable directory path | Raw External source catalog root. | Set when running `process_extcat`. | CLI override; rebuild not required |
| `EXTCAT_OUTPUT_DIRECTORY` | `const char*` | `LensingConfig::SOURCE_CAT_DEFAULT` | same | `--extcat-output` | Writable tile directory | Tile output and effective `process_main` External source catalog path. | Change for each catalog deployment. | No via CLI; yes after editing the header |
| `EXTCAT_FILENAME_TOKENS` | `std::vector<std::string>` | empty | empty | Repeatable `--extcat-contains` | Non-empty basename tokens; OR matching | Filters raw catalog files. | Change when filenames need filtering. | CLI override; rebuild not required |
| `EXTCAT_RECURSIVE` | `bool` | `true` | `true` | `--extcat-recursive` | Boolean | Recurse below the raw catalog root. | Disable for a flat directory only. | CLI override; rebuild not required |
| `EXTCAT_DELIMITER` | `const char*` | `"auto"` | same | `--extcat-delimiter` | `auto`, `whitespace`, `comma`, `tab` | Raw table delimiter mode. | Change when auto-detection is unsuitable. | CLI override; rebuild not required |
| `EXTCAT_HEADER_MODE` | `const char*` | `"auto"` | same | `--extcat-header` | `auto`, `present`, `absent` | Raw table header handling. | Change when auto-detection is unsuitable. | CLI override; rebuild not required |
| `EXTCAT_MALFORMED_POLICY` | `const char*` | `"fail"` | same | `--extcat-malformed` | `fail`, `skip` | Malformed-row policy. | Use `skip` only after accepting data loss. | CLI override; rebuild not required |
| `EXTCAT_EXISTING_POLICY` | `const char*` | `"fail"` | same | `--extcat-existing` | `fail`, `overwrite` | Existing-tile policy. | Select intentionally for reruns. | CLI override; rebuild not required |
| `EXTCAT_CHUNK_MIB` | `std::uint64_t` | `64` | `64` | `--extcat-chunk-mib` | Positive MiB value | MPI byte-range task size. | Tune for storage and rank count. | CLI override; rebuild not required |
| `EXTCAT_TOTAL_COLUMNS` | `std::size_t` | `18` | `18` | No | Positive canonical width | Pass-through External source catalog width and downstream layout basis. | Change only with coordinated schema consumers. | Yes |
| `EXTCAT_USE_EXPLICIT_COLUMNS` | `bool` | `false` | `false` | Set true by `--extcat-columns` | Boolean | Enables ordered column projection. | Use when raw tables contain extra/reordered fields. | CLI override; rebuild not required |
| `EXTCAT_INPUT_COLUMNS_ONE_BASED` | `std::vector<std::size_t>` | `{1, ..., 18}` | same | `--extcat-columns` | Non-empty positive one-based indices | Ordered raw columns emitted to tiles. | Change for another raw schema. | CLI override; rebuild not required |
| `EXTCAT_USE_EXPLICIT_COORDINATE_COLUMNS` | `bool` | `false` | `false` | Set true by any coordinate-column CLI | Boolean | Bypasses header coordinate discovery. | Use for known positional schemas. | CLI override; rebuild not required |
| `EXTCAT_RA_COLUMN_ONE_BASED` | `std::size_t` | `5` | `5` | `--extcat-ra-column` | Positive one-based index | Raw `ra` field position. | Change for another catalog schema. | CLI override; rebuild not required |
| `EXTCAT_DEC_COLUMN_ONE_BASED` | `std::size_t` | `6` | `6` | `--extcat-dec-column` | Positive one-based index | Raw `dec` field position. | Change for another catalog schema. | CLI override; rebuild not required |
| `EXTCAT_ZP_COLUMN_ONE_BASED` | `std::size_t` | `17` | `17` | `--extcat-zp-column` | Positive one-based index | Raw `zp` field position consumed downstream. | Change for another catalog schema. | CLI override; rebuild not required |

## `LensingConfig` (`config/LensingConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `pi` | `double` | `3.14159265358979323846` | same | No | Mathematical constant | Supplies angle conversions. | Derived parameter — do not edit directly. | Yes |
| `arc_convert` | `double` | `pi / 180` | same | No | Radians per degree | Converts degrees to radians. | Derived parameter — do not edit directly. | Yes |
| `npx` | `int` | `3000` | same | No | Positive pixels | Nominal CCD image width. | Change only for another detector contract. | Yes |
| `npy` | `int` | `5000` | same | No | Positive pixels | Nominal CCD image height. | Change only for another detector contract. | Yes |
| `ASTROMETRY_trivial` | `int` | `0` | `N/A — removed in Lite` | No | `0` Gaia, `1` identity | Selects Standard astrometry branch. | Debug or deliberately bypass Gaia only. | Yes |
| `AstroCatType` | `int` | `1` | same | No | `1` = legacy large Gaia tiles; `2` = one-degree Type-2 tiles | Selects the Stage-1 Gaia filename/read layout without changing the two-column row schema. | Set to `2` when `ASTROMETRY_CAT` points to `process_astrocat` output. | Yes |
| `PROCESS_stage` | `int` | `223092870` | same | No | Product of stage primes `2,3,5,7,11,13,17,19,23`; 23 requires 19 | Enables numerical stages by divisibility. | Change for staged/restart runs. | Yes |
| `include_FLAT` | `int` | `0` | `N/A — removed in Lite` | No | `0` off, `1` on | Enables Standard super-flat correction. | Enable only with valid flat files. | Yes |
| `include_Mask` | `int` | `2` | `N/A — removed in Lite` | No | `0` none, `1` legacy, `2` per-chip DQ, `3` both | Selects Standard mask branch; Lite is fixed to per-chip DQ. | Change when DQ masks are unavailable or mask mode changes. | Yes |
| `include_BGsub` | `int` | `1` | same | No | `0` off, `1` on | Subtracts the fitted Science-image background. | Change only for controlled preprocessing experiments. | Yes |
| `ASTROMETRY_CAT` | `std::string` | `/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted` | same | No | Readable Gaia tile directory | Gaia catalog location. | Change for every site/catalog deployment. | Yes |
| `SOURCE_CAT_DEFAULT` | `const char*` | `/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat` | same | `--extcat-output` | Readable/writable tile directory | Compiled External source catalog default. | Prefer CLI per deployment. | CLI override; rebuild not required |
| `FLAT_PATH` | `std::string` | `/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014` | `N/A — removed in Lite` | No | Readable flat FITS directory | Standard super-flat location. | Change when `include_FLAT=1`. | Yes |
| `PSF_PATH` | `std::string` | `"hahahaha"` | `N/A — removed in Lite` | No | Readable PSF image directory | Standard external-PSF location. | Replace before `ext_PSF=1`. | Yes |
| `ext_cat` | `int` | `1` | `N/A — removed in Lite` | No | `0` off, `1` on | Selects Standard External source catalog branch; Lite is fixed on. | Change only for a Standard no-catalog run. | Yes |
| `ext_PSF` | `int` | `0` | `N/A — removed in Lite` | No | `0` frame stars, `1` external PSF | Selects Standard PSF source; Lite uses frame stars. | Change only with valid external PSFs. | Yes |
| `CCD_split` | `int` | `2` | same | No | `1` whole chip, `2` amplifier split | Sets background/noise amplifier regions. | Change for another detector/readout model. | Yes |
| `nct` | `int` | `12` | same | No | Positive rectangle count | Number of background rectangles. | Tune only with background-model validation. | Yes |
| `ncx` | `int` | `3` | same | No | Positive x count | Background rectangles along x. | Keep consistent with `nct` geometry. | Yes |
| `psf_order` | `int` | `8` | same | No | Supported PSF polynomial selector | Exposure PSF polynomial order. | Adjust PSF modeling only. | Yes |
| `npo` | `int` | `64` | same | No | Positive sample count | Exposure PSF sample count. | Adjust PSF sampling only. | Yes |
| `npox` | `int` | `8` | same | No | Positive x count | Exposure PSF samples along x. | Keep consistent with `npo`. | Yes |
| `nstar_min` | `int` | `npo * 3 / 2` = `96` | same | No | Derived positive count | Minimum stars for exposure PSF fitting. | Derived parameter — do not edit directly. | Yes |
| `npl` | `int` | `10` | same | No | Non-negative coefficient-count offset | Local PSF coefficient count minus one. | Adjust local PSF model only. | Yes |
| `nplx` | `int` | `2` | same | No | Non-negative x degree | Local PSF polynomial x degree. | Keep consistent with `npl`. | Yes |
| `nstar_min_local` | `int` | `16` | same | No | Positive count | Minimum stars retained for a local fit. | Tune only with PSF fit validation. | Yes |
| `psf_exposure_min_candidates` | `int` | `60` | same | No | Positive count | Minimum exposure-wide PSF candidates. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_hist_bins` | `int` | `128` | same | No | Integer ≥ 3 | FWHM-locus histogram bins. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_locus_sigma` | `double` | `4.0` | same | No | Positive sigma multiplier | Exposure FWHM-locus window. | Tune Stage 5 selection. | Yes |
| `psf_fwhm_locus_min_samples` | `int` | `30` | same | No | Positive count | Minimum samples for a valid FWHM locus. | Tune sparse-exposure handling. | Yes |
| `PsfGroupingType` | `int` | `2` | same | No | `1` threshold graph, `2` mutual KNN | Selects PSF grouping topology. | Change for controlled algorithm comparison. | Yes |
| `psf_minchi_reference_fraction` | `double` | `1 / 3` | same | No | `(0, 1]` | Fraction of largest locus candidates eligible as references. | Tune Stage 5 threshold estimation. | Yes |
| `psf_minchi_reference_max_per_chip` | `int` | `5` | same | No | Positive count | Caps reference stars per chip. | Tune Stage 5 threshold estimation. | Yes |
| `psf_minchi_sigma_cut` | `double` | `4.0` | same | No | Positive sigma multiplier | Minimum-chi upper-tail rejection cut. | Tune Stage 5 selection. | Yes |
| `psf_knn_k` | `int` | `8` | same | No | Positive neighbor count | Neighbors in mutual-KNN grouping. | Change with grouping validation. | Yes |
| `psf_group_merge_ratio` | `double` | `0.30` | same | No | Non-negative group-size ratio | Secondary/main group merge threshold. | Tune Stage 5 grouping. | Yes |
| `psf_group_merge_min_gaia` | `int` | `2` | same | No | Non-negative match count | Gaia support required to merge a secondary group. | Tune Stage 5 grouping. | Yes |
| `psf_gaia_match_radius_pix` | `double` | `2.5` | same | No | Positive pixels | Gaia matching radius for PSF candidates. | Change for astrometric precision/pixel scale. | Yes |
| `psf_gaia_locus_min_matches` | `int` | `10` | same | No | Positive match count | Gaia matches required for locus support. | Tune sparse fields. | Yes |
| `psf_press_rejection_enabled` | `bool` | `true` | same | No | Boolean | Enables optional standardized-PRESS cleanup. | Disable for controlled fallback testing. | Yes |
| `psf_press_sigma_cut` | `double` | `4.0` | same | No | Positive sigma multiplier | PRESS outlier rejection cut. | Tune only with PSF residual validation. | Yes |
| `psf_press_max_removals` | `int` | `5` | same | No | Non-negative count | Maximum proposed PRESS removals per chip. | Tune only with PSF residual validation. | Yes |
| `psf_loo_min_denom` | `double` | `1.0e-6` | same | No | `(0, 1)` | Minimum analytic leave-one-out denominator. | Numerical guard; normally unchanged. | Yes |
| `step_psf` | `int` | `100` | `N/A — removed in Lite` | No | Positive pixels | Standard PSF-star spatial sampling step. | Change only for the corresponding Standard branch. | Yes |
| `deblending` | `int` | `1` | `N/A — removed in Lite` | No | `0` off, `1` on | Standard source deblending selector; Lite is fixed on. | Change only for controlled source tests. | Yes |
| `n_neighbor` | `int` | `5` | `N/A — removed in Lite` | No | Positive neighbor count | Standard deblending neighborhood size. | Change with deblending validation. | Yes |
| `PSF_type` | `int` | `1` | `N/A — removed in Lite` | No | `1` local polynomial, `2` hybrid | Standard PSF model selector; Lite is fixed to local. | Change only for Standard hybrid PSF. | Yes |
| `PSF_Ms` | `int` | `0` | `N/A — removed in Lite` | No | `0` local only, `1` PCA/multi-scale | Standard PCA branch selector; Lite has no PCA branch. | Change only with PCA inputs/resources. | Yes |
| `ns` | `int` | `64` | same | No | Positive even pixels | Science stamp and Fourier-grid side. | Change only with coordinated stamp/FFT validation. | Yes |
| `nsns` | `int` | `ns * ns` = `4096` | same | No | Derived pixels | Pixels in one science stamp. | Derived parameter — do not edit directly. | Yes |
| `chip_margin` | `int` | `8` | same | No | Non-negative pixels | Extra extraction margin near chip edges. | Change with stamp geometry. | Yes |
| `ns_2` | `int` | `ns / 2` = `32` | same | No | Derived pixels | Half stamp side. | Derived parameter — do not edit directly. | Yes |
| `nl_2` | `int` | `ns_2 + chip_margin` = `40` | same | No | Derived pixels | Half expanded extraction side. | Derived parameter — do not edit directly. | Yes |
| `nl` | `int` | `nl_2 * 2` = `80` | same | No | Derived pixels | Full expanded extraction side. | Derived parameter — do not edit directly. | Yes |
| `flag_thresh` | `int` | `3` | same | No | Non-negative extraction flag | Maximum accepted extraction flag. | Tune source quality selection. | Yes |
| `chip_edge_margin` | `int` | `chip_margin` = `8` | same | No | Derived pixels | Alias used by chip-edge checks. | Derived parameter — do not edit directly. | Yes |
| `dz_thresh` | `double` | `0.1` | same | No | Non-negative redshift difference | Redshift tolerance for catalog matching. | Tune matching for another catalog/error model. | Yes |
| `len_g` | `int` | `40` | same | No | Positive internal row width | Galaxy metadata capacity. | Internal layout; normally unchanged. | Yes |
| `len_s` | `int` | `15` | same | No | Positive internal row width | Star metadata capacity. | Internal layout; normally unchanged. | Yes |
| `n_user_max` | `int` | `200` | same | No | Positive count | Bright detections used for astrometric matching. | Tune astrometric pattern matching. | Yes |
| `ngal_max` | `int` | `4000` | same | No | Positive reservation hint | Initial galaxy-vector capacity; not a catalog limit. | Change only for allocation tuning. | Yes |
| `nstar_max` | `int` | `2000` | same | No | Positive reservation hint | Initial star-vector capacity; not a catalog limit. | Change only for allocation tuning. | Yes |
| `npara` | `int` | `25` | same | No | Fixed output width | Number of pipeline fields through exposure `chi2`. | Internal schema; coordinate all readers/writers. | Yes |
| `len_sam` | `int` | `50` | same | No | Positive internal row width | PSF sample metadata length. | Internal layout; normally unchanged. | Yes |
| `npd` | `int` | `33` | same | No | Positive coefficient count | PU astrometric distortion terms. | Change only with astrometric model code. | Yes |
| `blocksize` | `int` | `200` | same | No | Positive pixels | Target background block side. | Tune background modeling. | Yes |
| `bg_rough_grid_x` | `int` | `32` | same | No | Positive grid count | Rough background grid columns. | Tune background modeling. | Yes |
| `bg_rough_grid_y` | `int` | `32` | same | No | Positive grid count | Rough background grid rows. | Tune background modeling. | Yes |
| `bg_min_block_pixels` | `int` | `1000` | same | No | Positive pixel count | Minimum pixels in a background block. | Tune masked/small images. | Yes |
| `bg_min_clipped_pixels` | `int` | `200` | same | No | Positive pixel count | Minimum pixels after block clipping. | Tune masked/small images. | Yes |
| `bg_min_valid_frac` | `double` | `0.25` | same | No | Fraction in `(0, 1]` | Minimum valid fraction per background block. | Tune masking tolerance. | Yes |
| `bg_clip_low` | `double` | `4.0` | same | No | Positive sigma multiplier | Lower background clipping limit. | Tune background robustness. | Yes |
| `bg_clip_high` | `double` | `2.5` | same | No | Positive sigma multiplier | Upper background clipping limit. | Tune background robustness. | Yes |
| `bg_fit_clip_sigma` | `double` | `3.0` | same | No | Positive sigma multiplier | Background-plane fit clipping. | Tune background robustness. | Yes |
| `bg_fit_max_iter` | `int` | `4` | same | No | Non-negative iterations | Maximum background-plane clipping iterations. | Tune convergence/runtime. | Yes |
| `bg_min_fit_factor` | `int` | `3` | same | No | Positive samples-per-coefficient factor | Minimum plane-fit sample multiplier. | Tune fit stability. | Yes |
| `source_thresh` | `double` | `2.0` | same | No | Positive S/N threshold | Source detection threshold. | Adjust scientific source selection. | Yes |
| `core_thresh` | `double` | `4.0` | same | No | Positive threshold | Source-core detection threshold. | Adjust scientific source selection. | Yes |
| `flat_thresh` | `double` | `0.01` | `N/A — removed in Lite` | No | Positive flat value | Minimum accepted Standard flat-field value. | Change with flat calibration. | Yes |
| `NstampType` | `int` | `1` | same | No | `1` physical blank stamp, `2` local covariance power | Selects Stage-3 noise product. | Change for controlled noise-method runs. | Yes |
| `noise_sigma_ratio_min` | `double` | `0.80` | same | No | Positive lower ratio | Blank/source sigma lower gate. | Tune blank-stamp quality. | Yes |
| `noise_sigma_ratio_max` | `double` | `1.25` | same | No | Above minimum | Blank/source sigma upper gate. | Tune blank-stamp quality. | Yes |
| `noise_mad_ratio_min` | `double` | `0.70` | same | No | Positive lower ratio | Blank/source MAD lower gate. | Tune blank-stamp quality. | Yes |
| `noise_mad_ratio_max` | `double` | `1.30` | same | No | Above minimum | Blank/source MAD upper gate. | Tune blank-stamp quality. | Yes |
| `noise_tail_sigma` | `double` | `2.5` | same | No | Positive sigma | Tail-count threshold. | Tune blank-stamp quality. | Yes |
| `noise_max_tail_fraction` | `double` | `0.05` | same | No | Fraction `[0, 1]` | Maximum blank-stamp tail fraction. | Tune blank-stamp quality. | Yes |
| `noise_max_mask_fraction` | `double` | `0.02` | same | No | Fraction `[0, 1]` | Maximum blank-stamp masked fraction. | Tune blank-stamp quality. | Yes |
| `noise_region_size` | `int` | `192` | same | No | Positive even pixels; greater than inner size | Outer local-noise square side. | Change with covariance geometry. | Yes |
| `noise_inner_size` | `int` | `96` | same | No | Even pixels; at least `nl` | Central covariance exclusion side. | Change with source/stamp geometry. | Yes |
| `noise_plane_min_valid_fraction` | `double` | `0.30` | same | No | Fraction `(0, 1]` | Minimum valid plane-fit shell fraction. | Tune masking tolerance. | Yes |
| `noise_cov_padding_factor` | `double` | `2.0` | same | No | Positive; padded side must support linear autocorrelation | Covariance FFT padding multiplier. | Change only with FFT validation. | Yes |
| `noise_cov_fft_size` | `int` | `384` | same | No | Derived padded side | Covariance FFT side. | Derived parameter — do not edit directly. | Yes |
| `noise_cov_max_lag` | `int` | `8` | same | No | `0 <= lag < noise_region_size` | Maximum retained signed covariance lag. | Tune covariance model. | Yes |
| `noise_cov_min_valid_pixels` | `int` | `4096` | same | No | Positive count | Minimum covariance-mask pixels. | Tune masking tolerance. | Yes |
| `noise_cov_min_pair_fraction` | `double` | `0.50` | same | No | Fraction `(0, 1]` | Minimum lag pair-count fraction. | Tune covariance reliability. | Yes |
| `noise_cov_sigma_ratio_min` | `double` | `0.80` | same | No | Positive lower ratio | Covariance/source sigma lower gate. | Tune covariance quality. | Yes |
| `noise_cov_sigma_ratio_max` | `double` | `1.25` | same | No | Above minimum | Covariance/source sigma upper gate. | Tune covariance quality. | Yes |
| `noise_cov_max_negative_fraction` | `double` | `0.25` | same | No | Fraction `[0, 1]` | Maximum negative power fraction. | Tune covariance quality. | Yes |
| `noise_cov_imag_tolerance` | `double` | `1.0e-10` | same | No | Non-negative numerical tolerance | Maximum imaginary FFT residual. | Numerical guard; normally unchanged. | Yes |
| `sig_blocksize` | `int` | `200` | same | No | Positive pixels | Mode-bar noise-estimator block side. | Tune only with noise-estimator calibration. | Yes |
| `sig_block_max` | `int` | `sig_blocksize²` = `40000` | same | No | Derived pixel count | Maximum pixels per noise block. | Derived parameter — do not edit directly. | Yes |
| `sig_max_blocks` | `int` | `2048` | same | No | Positive count | Maximum sampled noise blocks. | Tune memory/runtime only. | Yes |
| `sig_min_block_pixels` | `int` | `1000` | same | No | Positive count | Minimum pixels in one noise block. | Tune sparse/masked data. | Yes |
| `sig_min_block_triples` | `int` | `1000` | same | No | Positive count | Minimum valid triples per block. | Tune sparse/masked data. | Yes |
| `sig_min_blocks` | `int` | `4` | same | No | Positive count | Minimum blocks for a plane fit. | Tune sparse data only. | Yes |
| `sig_hist_nbin` | `int` | `256` | same | No | Positive bin count | Mode-finding histogram bins. | Tune estimator resolution. | Yes |
| `sig_hist_range` | `double` | `6.0` | same | No | Positive sigma range | Mode histogram range. | Tune estimator robustness. | Yes |
| `sig_min_mode_count` | `int` | `500` | same | No | Positive count | Minimum samples defining the mode. | Tune sparse data only. | Yes |
| `sig_min_lower_count` | `int` | `1000` | same | No | Positive count | Minimum lower-side width samples. | Tune sparse data only. | Yes |
| `sig_lower_quantile` | `double` | `0.3173105` | same | No | Quantile in `(0, 1)` | Lower-side width quantile. | Calibration constant; normally unchanged. | Yes |
| `sig_clip_k` | `double` | `3.0` | same | No | Positive sigma multiplier | Symmetric clipping threshold. | Tune estimator robustness. | Yes |
| `sig_rdil` | `int` | `2` | same | No | Positive pixel stride | Pixel stride used by the estimator. | Tune sampling/runtime only. | Yes |
| `sig_clip_niter` | `int` | `2` | same | No | Non-negative iterations | Number of clipping iterations. | Tune convergence/runtime. | Yes |
| `sig_min_fit_triples` | `int` | `1000` | same | No | Positive count | Minimum triples in final fit. | Tune sparse data only. | Yes |
| `sig_min_fit_frac` | `double` | `0.20` | same | No | Fraction `(0, 1]` | Minimum retained fit fraction. | Tune robustness only. | Yes |
| `sig_median_ratio` | `double` | `1.2678405` | same | No | Positive calibration factor | Median-to-sigma conversion. | Calibration constant; normally unchanged. | Yes |
| `sig_plane_min` | `double` | `1.0e-8` | same | No | Positive floor | Minimum noise-plane value. | Numerical guard; normally unchanged. | Yes |
| `sig_max_plane_ratio` | `double` | `4.0` | same | No | Ratio ≥ 1 | Maximum noise-plane variation. | Tune rejection only with validation. | Yes |
| `sig_pivot_min` | `double` | `1.0e-8` | same | No | Positive floor | Minimum linear-solve pivot. | Numerical guard; normally unchanged. | Yes |
| `sig_scale_s1` | `double` | `0.673475` | same | No | Positive calibration candidate | Stage-1 noise calibration candidate. | Calibration experiments only. | Yes |
| `sig_scale_s2` | `double` | `1.027786` | same | No | Positive calibration value | Stage-2 noise calibration. | Calibration experiments only. | Yes |
| `sig_scale` | `double` | `sig_scale_s2` | same | No | Derived active selector | Active noise calibration scale. | Derived parameter — select a calibrated source value. | Yes |
| `area_max` | `int` | `ns²` = `4096` | same | No | Derived pixels | Maximum connected source area. | Derived parameter — do not edit directly. | Yes |
| `area_thresh` | `int` | `6` | same | No | Positive pixels | Minimum connected source area. | Adjust scientific source selection. | Yes |
| `gal_smooth` | `int` | `0` | same | No | Supported smoothing selector | Galaxy-stamp smoothing type. | Change for controlled processing tests. | Yes |
| `star_smooth` | `int` | `2` | same | No | Supported smoothing selector | Star-stamp smoothing type. | Change for controlled PSF tests. | Yes |
| `SNR_PSF` | `double` | `100.0` | same | No | Positive S/N threshold | Minimum PSF-star S/N. | Adjust PSF star selection. | Yes |
| `saturation_thresh` | `double` | `25000.0` | same | No | Detector-count threshold | Saturated-pixel threshold. | Change for detector/gain regime. | Yes |
| `pixel_size` | `double` | `0.2628` | same | No | Positive arcsec/pixel | DECam pixel scale. | Change only for another detector. | Yes |
| `iid` | `int` | `0` | same | No | Derived zero-based index | PSF polynomial `chi2` output column. | Derived schema — do not edit directly. | Yes |
| `ipixx` | `int` | `1` | same | No | Derived zero-based index | Source-center x column. | Derived schema — do not edit directly. | Yes |
| `ipixy` | `int` | `2` | same | No | Derived zero-based index | Source-center y column. | Derived schema — do not edit directly. | Yes |
| `isig` | `int` | `3` | same | No | Derived zero-based index | Local noise sigma column. | Derived schema — do not edit directly. | Yes |
| `istar` | `int` | `4` | same | No | Derived zero-based index | Available PSF-star count column. | Derived schema — do not edit directly. | Yes |
| `ipeak` | `int` | `4` | same | No | Derived legacy alias | Historical peak alias. | Derived schema — do not edit directly. | Yes |
| `i_imax` | `int` | `5` | same | No | Derived zero-based index | Peak x column. | Derived schema — do not edit directly. | Yes |
| `i_jmax` | `int` | `6` | same | No | Derived zero-based index | Peak y column. | Derived schema — do not edit directly. | Yes |
| `ih_flux` | `int` | `7` | same | No | Derived zero-based index | Half-light flux column. | Derived schema — do not edit directly. | Yes |
| `ih_area` | `int` | `8` | same | No | Derived zero-based index | Source area column. | Derived schema — do not edit directly. | Yes |
| `iflag` | `int` | `9` | same | No | Derived zero-based index | Quality flag column. | Derived schema — do not edit directly. | Yes |
| `iPSF` | `int` | `10` | same | No | Derived zero-based index | Local PSF size column. | Derived schema — do not edit directly. | Yes |
| `iSNR_F` | `int` | `11` | same | No | Derived zero-based index | Fourier S/N column. | Derived schema — do not edit directly. | Yes |
| `ira` | `int` | `12` | same | No | Derived zero-based index | Source RA column. | Derived schema — do not edit directly. | Yes |
| `idec` | `int` | `13` | same | No | Derived zero-based index | Source Dec column. | Derived schema — do not edit directly. | Yes |
| `igf1` | `int` | `14` | same | No | Derived zero-based index | Field-distortion `g1` column. | Derived schema — do not edit directly. | Yes |
| `igf2` | `int` | `15` | same | No | Derived zero-based index | Field-distortion `g2` column. | Derived schema — do not edit directly. | Yes |
| `ig1` | `int` | `16` | same | No | Derived zero-based index | Fourier_Quad `g1` column. | Derived schema — do not edit directly. | Yes |
| `ig2` | `int` | `17` | same | No | Derived zero-based index | Fourier_Quad `g2` column. | Derived schema — do not edit directly. | Yes |
| `ide` | `int` | `18` | same | No | Derived zero-based index | Shear response column. | Derived schema — do not edit directly. | Yes |
| `ih1` | `int` | `19` | same | No | Derived zero-based index | Higher-order `h1` column. | Derived schema — do not edit directly. | Yes |
| `ih2` | `int` | `20` | same | No | Derived zero-based index | Higher-order `h2` column. | Derived schema — do not edit directly. | Yes |
| `icos2` | `int` | `21` | same | No | Derived zero-based index | Spin-2 cosine column. | Derived schema — do not edit directly. | Yes |
| `isin2` | `int` | `22` | same | No | Derived zero-based index | Spin-2 sine column. | Derived schema — do not edit directly. | Yes |
| `iparity` | `int` | `23` | same | No | Derived zero-based index | WCS parity column. | Derived schema — do not edit directly. | Yes |
| `ichi2` | `int` | `24` | same | No | Derived zero-based index | Exposure `chi2` column. | Derived schema — do not edit directly. | Yes |
| `NMAX_EXPO` | `int` | `25000` | same | No | Positive exposure count | Maximum exposures accepted per run. | Change only for larger validated workloads. | Yes |
| `NMAX_CHIP` | `int` | `62` | same | No | Positive chip count | Maximum CCDs per exposure. | Change only for another camera/layout. | Yes |
| `g1_c` | `double` | `-0.001` | same | No | Additive calibration | Field-distortion `g1` correction. | Recalibrate for another dataset/band. | Yes |
| `g2_c` | `double` | `-0.0003` | same | No | Additive calibration | Field-distortion `g2` correction. | Recalibrate for another dataset/band. | Yes |
| `chi2_thresh` | `double` | `0.01` | same | No | Non-negative threshold | Maximum exposure PSF `chi2`. | Adjust scientific quality selection. | Yes |
| `chipnx` | `int` | `2046` | same | No | Positive pixels | Science CCD width used for PSF coordinates. | Change only for another detector. | Yes |
| `chipny` | `int` | `4094` | same | No | Positive pixels | Science CCD height used for PSF coordinates. | Change only for another detector. | Yes |
| `rescale_size` | `double` | `1.2` | `N/A — removed in Lite` | No | Positive target size | Standard PCA PSF residual rescaling size. | Change only for `PSF_Ms=1`. | Yes |
| `procs_pn` | `int` | `40` | `N/A — removed in Lite` | No | Positive rank count | MPI ranks per PCA scheduling group. | Tune Standard PCA scheduling. | Yes |
| `work_pn` | `int` | `10` | `N/A — removed in Lite` | No | Positive worker count | Concurrent PCA workers per group. | Tune Standard PCA scheduling. | Yes |
| `nblocks` | `int` | `2` | `N/A — removed in Lite` | No | Positive blocks per CCD axis | Standard PCA spatial grid. | Tune Standard PCA modeling. | Yes |
| `n_pcs` | `int` | `100` | `N/A — removed in Lite` | No | Positive component count | Maximum PCA components. | Tune Standard PCA modeling. | Yes |
| `npp6th` | `int` | `28` | `N/A — removed in Lite` | No | Polynomial term count | Sixth-degree 2D PCA surface terms. | Derived by polynomial basis; normally unchanged. | Yes |
| `pca_negative_eigenvalue_threshold` | `double` | `-1.0e-5` | `N/A — removed in Lite` | No | Non-positive tolerance | Invalid PCA eigenvalue cutoff. | Numerical guard; normally unchanged. | Yes |
| `nmax_star_pchip` | `int` | `1000000` | `N/A — removed in Lite` | No | Positive reservation capacity | Legacy PCA star capacity. | Change only for Standard PCA memory planning. | Yes |

## `ProcessRearrConfig` (`config/ProcessRearrConfig.hpp` and `config/pathconfig.hpp`)

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `ichi2` | `std::size_t` | `25` | `25` | No | Derived from `LensingConfig::ichi2 + 1` | Number of pipeline fields appended after the identity columns. | Derived parameter — do not edit directly. | Yes |
| `EXPO_COLUMN_COUNT` | `std::size_t` | `1` | `1` | No | Fixed schema count | Number of original-exposure fields. | Derived schema — do not edit directly. | Yes |
| `CCD_COLUMN_COUNT` | `std::size_t` | `1` | `1` | No | Fixed schema count | Number of CCD-number fields. | Derived schema — do not edit directly. | Yes |
| `ALL_CAT_TOTAL_COLUMNS` | `std::size_t` | `45` | `45` | No | `18 + 1 + 1 + 25` | Default complete `_all.cat` row width. | Derived parameter — do not edit directly. | Yes |
| `SKY_GRID_DEGREES` | `double` | `0.1` | same | No | Positive degrees | Full-sky RA/Dec grid width. | Change for another spatial partition resolution. | Yes |
| `RA_BIN_COUNT` | `int` | `3600` | same | No | Positive full-sky bin count | RA grid dimension. | Keep consistent with grid width. | Yes |
| `DEC_BIN_COUNT` | `int` | `1800` | same | No | Positive full-sky bin count | Dec grid dimension. | Keep consistent with grid width. | Yes |
| `SKY_TILE_COUNT` | `std::size_t` | `6480000` | same | No | `RA_BIN_COUNT * DEC_BIN_COUNT` | Total full-sky grid cells. | Derived parameter — do not edit directly. | Yes |
| `TARGET_SUBCAT_ROWS` | `std::uint64_t` | `500000` | same | No | Positive row count | Target rows per weighted k-d partition. | Tune output file size and memory. | Yes |
| `SKIP_DIRECTORY_NAME` | `std::string_view` | `"Large_Field"` | same | No | Directory name | Directory excluded from catalog scans. | Change only when layout conventions change. | Yes |
| `SUBCAT_PREFIX` | `std::string_view` | `"subcat_"` | same | No | Filename prefix | Spatial partition filename prefix. | Change for another naming convention. | Yes |
| `SUBCAT_EXTENSION` | `std::string_view` | `".cat"` | same | No | Filename extension | Spatial partition filename extension. | Change for another naming convention. | Yes |
| `SUBCAT_ID_WIDTH` | `int` | `6` | same | No | Positive digit count | Minimum zero-padded partition-ID width. | Change for another naming convention. | Yes |
| `SUMMARY_FILENAME` | `std::string_view` | `"catalog_summary.txt"` | same | No | Filename | Rearrangement summary name. | Change for another naming convention. | Yes |
| `OUTPUT_PRECISION` | `int` | `10` | same | No | Positive significant digits | Precision of output catalog values. | Change only after precision/size review. | Yes |
| `SUMMARY_PRECISION` | `int` | `4` | same | No | Non-negative decimal places | Precision of summary bounds. | Change for reporting needs. | Yes |
| `SKIP_MISSING_CATALOGS` | `bool` | `true` | same | No | Boolean | Continues past missing `_all.cat` inputs. | Set false for strict completeness checks. | Yes |
| `SKIP_MALFORMED_ROWS` | `bool` | `true` | same | No | Boolean | Continues past malformed rows. | Set false for strict schema checks. | Yes |
| `externalCatalogColumns(options)` | function result | `18` without projection | same | `--extcat-columns` changes result | Projection length or `EXTCAT_TOTAL_COLUMNS` | Runtime-effective external prefix width. | Derived parameter — do not edit directly. | No |
| `allCatalogColumns(options)` | function result | `45` without projection | same | `--extcat-columns` changes result | External width + `2` + `25` | Runtime-effective complete row width. | Derived parameter — do not edit directly. | No |

## `OutputLayout` (`config/pathconfig.hpp`)

These arrays are identical in Standard and Lite and have no CLI override. They
contain relative directory names, not deployment roots; the runtime dataset root
is prepended by the existing path helpers.

| Parameter | Type | Standard / Lite compiled value | Function | Rebuild after change |
|---|---|---|---|---|
| `NON_CHIP_BASE_DIRECTORIES` | `std::array<const char*, 14>` | `science`, `dqmask`, `stamps`, `result`, `stamps/dat_StarInfo`, `stamps/fits_StarP`, `stamps/fits_PsfSrc`, `stamps/dat_ExpoInfo`, `stamps/dat_StarComp`, `stamps/dat_Rescale`, `stamps/dat_Pcs`, `stamps/dat_StarCompV2`, `astrometry/Head`, `astrometry/dat_Chk` | Complete fixed base-directory contract created without a chip suffix. | Yes |
| `CHIP_PRODUCT_DIRECTORIES` | `std::array<const char*, 16>` | `stamps/Norm`, `stamps/cat_Orig`, `stamps/dat_StarCanInfo`, `stamps/fits_StarCan`, `stamps/fits_StarCanN`, `stamps/fits_StarCanP`, `stamps/dat_SrcInfo`, `stamps/fits_Src`, `stamps/fits_Noise`, `stamps/fits_SrcP`, `stamps/dat_PsfFit`, `stamps/fits_PsfLocal`, `stamps/dat_Shear`, `stamps/dat_StarXY`, `stamps/fits_PsfResi`, `astrometry/dat_Astro` | Complete fixed per-chip product-directory contract. | Yes |

`include/general/OutputLayout.hpp` now contains only the functions that derive
exposure and chip paths from these centralized arrays.

## `config/FDConfig.hpp`

| Parameter | Type | Standard default | Lite default | CLI override | Legal values / meaning | Function | When to change | Rebuild after change |
|---|---|---|---|---|---|---|---|---|
| `FD_STATIC_MODE` | `StaticMode` | `PDF_SIGMA` | same | No | `PDF_SIGMA`, `PDF_JACK`, `SWSE_JACK` | Selects FD estimator/uncertainty mode. | Change for another validated statistical method. | Yes |
| `FD_USE_PDF_STATIS` | `bool` | `true` | same | No | Derived from `FD_STATIC_MODE` | Enables PDF sign-test path. | Derived parameter — do not edit directly. | Yes |
| `FD_USE_JACKKNIFE` | `bool` | `false` | same | No | Derived from `FD_STATIC_MODE` | Enables jackknife uncertainty. | Derived parameter — do not edit directly. | Yes |
| `FD_USE_SWSE_DATA` | `bool` | `false` | same | No | Derived from `FD_STATIC_MODE` | Enables SWSE data model. | Derived parameter — do not edit directly. | Yes |
| `FD_PER_EXPOSURE_STAR_BAR` | `bool` | `false` | same | No | `false` global N-to-1, `true` per-exposure N-to-N | Selects star-bar fitting scope. | Change for another validated FD mode. | Yes |
| `nmax_per_core` | `int` | `20000000` | same | No | Positive source capacity | Maximum sources reserved per MPI rank. | Tune memory/capacity for workload. | Yes |
| `fd_num` | `int` | `21` | same | No | Positive bin count | Spatial bins by field distortion. | Tune FD spatial resolution. | Yes |
| `PDF_BINS` | `int` | `4` | same | No | Positive inner-bin count | Equal-probability PDF bins. | Tune PDF estimator resolution. | Yes |
| `gf_lim` | `float` | `0.0015` | same | No | Positive distortion half-range | FD spatial bin range `±gf_lim`. | Change for another distortion range. | Yes |
| `NMAX` | `int` | `200` | same | No | Positive sampling count | Fine-grid sampling points. | Tune fit resolution/runtime. | Yes |
| `MAX_DUP` | `int` | `5` | same | No | Positive duplicate limit | Maximum duplicate measurements. | Change for another catalog duplication policy. | Yes |
| `N_jack` | `int` | `50` | same | No | Positive region count | Jackknife regions. | Tune covariance resolution. | Yes |
| `nmax_total` | `int` | `1000000` | same | No | Positive source count | Maximum sources used by K-means. | Tune memory/runtime. | Yes |
| `Km_iter` | `int` | `100` | same | No | Positive iterations | K-means iteration count. | Tune convergence/runtime. | Yes |
| `snrfcut` | `float` | `4.0` | same | No | Non-negative S/N | Minimum Fourier S/N cut. | Adjust scientific selection. | Yes |
| `snrlow` | `float` | `0.0` | same | No | `0` disables or lower S/N bound | Optional lower source-S/N bound. | Enable for a bounded S/N sample. | Yes |
| `snrhigh` | `float` | `0.0` | same | No | `0` disables or upper S/N bound | Optional upper source-S/N bound. | Enable for a bounded S/N sample. | Yes |
| `starcut` | `float` | `20.0` | same | No | Size threshold | Point-source size cut. | Recalibrate stellar selection. | Yes |
| `chi2_thresh` | `float` | `0.01` | same | No | Non-negative threshold | Maximum exposure `chi2`. | Adjust scientific quality selection. | Yes |
| `flagcut` | `float` | `0.0` | same | No | Maximum accepted flag | Source quality flag cut. | Adjust scientific quality selection. | Yes |
| `imaxcut` | `float` | `64.0` | same | No | Pixel-coordinate upper cut | Maximum peak x. | Change with stamp geometry. | Yes |
| `jmaxcut` | `float` | `64.0` | same | No | Pixel-coordinate upper cut | Maximum peak y. | Change with stamp geometry. | Yes |
| `zplow` | `float` | `0.0` | same | No | Lower `zp` bound | Minimum catalog `zp`. | Adjust redshift/sample selection. | Yes |
| `zphigh` | `float` | `3.0` | same | No | Upper `zp` bound above `zplow` | Maximum catalog `zp`. | Adjust redshift/sample selection. | Yes |
| `r_half_thresh` | `float` | `0.0` | same | No | `0` disables or size threshold | Optional half-light-radius threshold. | Enable for a size-selected sample. | Yes |
| `star_bar_mltp` | `float` | `3.0` | same | No | Positive sigma multiplier | Star-bar rejection multiplier. | Recalibrate stellar selection. | Yes |
| `psf_chi2_mltp` | `float` | `3.0` | same | No | Positive sigma multiplier | PSF `chi2` rejection multiplier. | Recalibrate PSF quality selection. | Yes |
| `ft_cut` | `float` | `-1.0` | same | No | Negative disables | `FLAGS_FT` cut. | Change for another External source catalog schema/selection. | Yes |
| `fg_cut` | `float` | `-10.0` | same | No | Catalog-specific cut | `FLAGS_FG` cut. | Change for another catalog selection. | Yes |
| `gold_cut` | `float` | `-10.0` | same | No | Catalog-specific cut | `FLAGS_GOLD` cut. | Change for another catalog selection. | Yes |
| `ext_cut` | `float` | `-4.0` | same | No | Catalog-specific cut | `EXT_MASH` cut. | Change for another catalog selection. | Yes |
| `n_size_bins` | `int` | `100` | same | No | Positive bin count | Stellar-size histogram bins. | Tune stellar-locus resolution. | Yes |
| `n_mag_bins` | `int` | `20` | same | No | Positive bin count | Magnitude histogram bins. | Tune stellar-locus resolution. | Yes |
| `size_min` | `float` | `-2.0` | same | No | Lower size bound | Stellar-size histogram minimum. | Recalibrate stellar locus. | Yes |
| `size_max` | `float` | `2.0` | same | No | Above `size_min` | Stellar-size histogram maximum. | Recalibrate stellar locus. | Yes |
| `mag_min_val` | `float` | `10.0` | same | No | Lower magnitude bound | Magnitude histogram minimum. | Change for another band/depth. | Yes |
| `mag_max_val` | `float` | `30.0` | same | No | Above `mag_min_val` | Magnitude histogram maximum. | Change for another band/depth. | Yes |
| `min_bin_count` | `int` | `100` | same | No | Positive sample count | Minimum samples in a usable bin. | Tune sparse samples. | Yes |
| `peak_match_tol` | `float` | `0.05` | same | No | Positive size tolerance | Stellar-peak matching tolerance. | Recalibrate stellar locus. | Yes |
| `min_concentration` | `float` | `0.6` | same | No | Fraction `[0, 1]` | Minimum stellar-locus concentration. | Recalibrate stellar locus. | Yes |
| `star_phy_min` | `float` | `-0.5` | same | No | Lower physical-size bound | Physical stellar-size minimum. | Recalibrate stellar locus. | Yes |
| `star_phy_max` | `float` | `0.2` | same | No | Above `star_phy_min` | Physical stellar-size maximum. | Recalibrate stellar locus. | Yes |
| `stage1_snr` | `float` | `40.0` | same | No | Non-negative S/N | Histogram-accumulation S/N. | Tune per-exposure star cuts. | Yes |
| `stage2_snr` | `float` | `0.0` | same | No | `0` disables or S/N threshold | Second-pass S/N threshold. | Enable only for a validated second pass. | Yes |
| `init_win_active` | `float` | `0.1` | same | No | Positive window | Active exposure initial size window. | Tune per-exposure star cuts. | Yes |
| `init_win_fallback` | `float` | `0.15` | same | No | Positive window | Fallback exposure initial size window. | Tune fallback behavior. | Yes |
| `default_s_init` | `float` | `0.5` | same | No | Size coordinate | Default initial stellar-size center. | Recalibrate stellar locus. | Yes |
| `clip_nsigma` | `float` | `3.0` | same | No | Positive sigma multiplier | Iterative clipping threshold. | Tune per-exposure robustness. | Yes |
| `min_clip_limit` | `float` | `0.015` | same | No | Positive half-width | Minimum clipping half-width. | Tune per-exposure robustness. | Yes |
| `default_s_std` | `float` | `0.05` | same | No | Positive scatter | Default stellar-size scatter. | Recalibrate stellar locus. | Yes |
| `fallback_scut_default` | `float` | `0.6` | same | No | Size threshold | Fallback stellar-size cut. | Recalibrate fallback behavior. | Yes |
| `col_flags_ft` | `int` | `0` | same | No | Derived zero-based index | External `FLAGS_FT` column. | Derived schema — do not edit directly. | Yes |
| `col_flags_fg` | `int` | `1` | same | No | Derived zero-based index | External `FLAGS_FG` column. | Derived schema — do not edit directly. | Yes |
| `col_flags_gold` | `int` | `2` | same | No | Derived zero-based index | External `FLAGS_GOLD` column. | Derived schema — do not edit directly. | Yes |
| `col_ext_mash` | `int` | `3` | same | No | Derived zero-based index | External `EXT_MASH` column. | Derived schema — do not edit directly. | Yes |
| `col_cra` | `int` | `4` | same | No | Derived zero-based index | External `ra` column. | Derived schema — do not edit directly. | Yes |
| `col_cdec` | `int` | `5` | same | No | Derived zero-based index | External `dec` column. | Derived schema — do not edit directly. | Yes |
| `col_mag_g` | `int` | `6` | same | No | Derived zero-based index | External g-band magnitude column. | Derived schema — do not edit directly. | Yes |
| `col_mag_r` | `int` | `8` | same | No | Derived zero-based index | External r-band magnitude column. | Derived schema — do not edit directly. | Yes |
| `col_mag_i` | `int` | `10` | same | No | Derived zero-based index | External i-band magnitude column. | Derived schema — do not edit directly. | Yes |
| `col_mag_z` | `int` | `12` | same | No | Derived zero-based index | External z-band magnitude column. | Derived schema — do not edit directly. | Yes |
| `col_mag_y` | `int` | `14` | same | No | Derived zero-based index | External y-band magnitude column. | Derived schema — do not edit directly. | Yes |
| `col_zp` | `int` | `16` | same | No | Derived zero-based index | External `zp` column. | Derived schema — do not edit directly. | Yes |
| `external_num_cols` | `int` | `18` | `18` | No | Optional in Standard, fixed in Lite | Effective external prefix width. | Derived schema — do not edit directly. | Yes |
| `col_expo` | `int` | `18` | same | No | Derived absolute index | Original 1-based exposure number. | Derived schema — do not edit directly. | Yes |
| `col_ccd` | `int` | `19` | same | No | Immediately follows `col_expo` | CCD-number column. | Derived schema — do not edit directly. | Yes |
| `source_col_offset` | `int` | `20` | same | No | External width + exposure + CCD | Prefix width before pipeline fields. | Derived parameter — do not edit directly. | Yes |
| `col_polychi2` | `int` | `20` | same | No | Derived absolute index | PSF fit `chi2` column. | Derived schema — do not edit directly. | Yes |
| `col_pixx` | `int` | `21` | same | No | Derived absolute index | Source x column. | Derived schema — do not edit directly. | Yes |
| `col_pixy` | `int` | `22` | same | No | Derived absolute index | Source y column. | Derived schema — do not edit directly. | Yes |
| `col_sig` | `int` | `23` | same | No | Derived absolute index | Noise sigma column. | Derived schema — do not edit directly. | Yes |
| `col_star` | `int` | `24` | same | No | Derived absolute index | PSF-star count column. | Derived schema — do not edit directly. | Yes |
| `col_peak` | `int` | `24` | same | No | Derived legacy alias | Historical peak column alias. | Derived schema — do not edit directly. | Yes |
| `col_imax` | `int` | `25` | same | No | Derived absolute index | Peak x column. | Derived schema — do not edit directly. | Yes |
| `col_jmax` | `int` | `26` | same | No | Derived absolute index | Peak y column. | Derived schema — do not edit directly. | Yes |
| `col_h_flux` | `int` | `27` | same | No | Derived absolute index | Half-light flux column. | Derived schema — do not edit directly. | Yes |
| `col_h_area` | `int` | `28` | same | No | Derived absolute index | Source area column. | Derived schema — do not edit directly. | Yes |
| `col_flag` | `int` | `29` | same | No | Derived absolute index | Quality flag column. | Derived schema — do not edit directly. | Yes |
| `col_PSF` | `int` | `30` | same | No | Derived absolute index | Local PSF size column. | Derived schema — do not edit directly. | Yes |
| `col_SNR_F` | `int` | `31` | same | No | Derived absolute index | Fourier S/N column. | Derived schema — do not edit directly. | Yes |
| `col_ra` | `int` | `32` | same | No | Derived absolute index | Source RA column. | Derived schema — do not edit directly. | Yes |
| `col_dec` | `int` | `33` | same | No | Derived absolute index | Source Dec column. | Derived schema — do not edit directly. | Yes |
| `col_gf1` | `int` | `34` | same | No | Derived absolute index | Field-distortion `g1` column. | Derived schema — do not edit directly. | Yes |
| `col_gf2` | `int` | `35` | same | No | Derived absolute index | Field-distortion `g2` column. | Derived schema — do not edit directly. | Yes |
| `col_g1` | `int` | `36` | same | No | Derived absolute index | Fourier_Quad `g1` column. | Derived schema — do not edit directly. | Yes |
| `col_g2` | `int` | `37` | same | No | Derived absolute index | Fourier_Quad `g2` column. | Derived schema — do not edit directly. | Yes |
| `col_de` | `int` | `38` | same | No | Derived absolute index | Shear response column. | Derived schema — do not edit directly. | Yes |
| `col_h1` | `int` | `39` | same | No | Derived absolute index | Higher-order `h1` column. | Derived schema — do not edit directly. | Yes |
| `col_h2` | `int` | `40` | same | No | Derived absolute index | Higher-order `h2` column. | Derived schema — do not edit directly. | Yes |
| `col_cos2` | `int` | `41` | same | No | Derived absolute index | Spin-2 cosine column. | Derived schema — do not edit directly. | Yes |
| `col_sin2` | `int` | `42` | same | No | Derived absolute index | Spin-2 sine column. | Derived schema — do not edit directly. | Yes |
| `col_parity` | `int` | `43` | same | No | Derived absolute index | WCS parity column. | Derived schema — do not edit directly. | Yes |
| `col_chi2` | `int` | `44` | same | No | Derived absolute index | Exposure `chi2` column. | Derived schema — do not edit directly. | Yes |
| `ICHI2` | `int` | `45` | same | No | Derived complete row width | Total default catalog columns. | Derived parameter — do not edit directly. | Yes |
| `bad_ccds` | `int[]` | `{2, 31, 53, 61}` | same | No | Detector-specific CCD IDs | DES CCDs excluded from FD analysis. | Change for another detector/quality list. | Yes |
| `n_bad_ccds` | `int` | `4` | same | No | Derived array length | Number of excluded CCDs. | Derived parameter — keep synchronized with `bad_ccds`. | Yes |
| `chip_xmin` | `int` | `50` | same | No | Pixel lower bound | Minimum accepted chip x. | Change for detector/edge-mask policy. | Yes |
| `chip_xmax` | `int` | `1990` | same | No | Pixel upper bound | Maximum accepted chip x. | Change for detector/edge-mask policy. | Yes |
| `chip_ymin` | `int` | `100` | same | No | Pixel lower bound | Minimum accepted chip y. | Change for detector/edge-mask policy. | Yes |
| `chip_ymax` | `int` | `3990` | same | No | Pixel upper bound | Maximum accepted chip y. | Change for detector/edge-mask policy. | Yes |

The default complete row is 18 External source catalog fields, one original
exposure number, one CCD number, and 25 pipeline fields: 45 columns total.
Explicit projection changes the
external prefix width used by `process_rearr`; the current FD column constants
remain tied to the default 18-field schema, so schema changes require coordinated
review of `ExtCatConfig`, `ExternalCatalogReader`, `ProcessRearrConfig`, and the FD
reader before rebuilding.
