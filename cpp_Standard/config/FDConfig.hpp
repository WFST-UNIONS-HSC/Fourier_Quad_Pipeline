#ifndef FD_CONFIG_HPP
#define FD_CONFIG_HPP

// ==========================================
// FDConfig — C++ equivalent of the Fortran para.inc
// Method: All FD-test parameters fixed to one set (DES defaults). The two
//         feature switches below select the sigma-estimation method and the
//         star-bar fitting mode at compile time.
// ==========================================

#include "LensingConfig.hpp"
#include "ExtCatConfig.hpp"

namespace FDConfig {

// ==================== Feature Switches ====================
// Statistical mode: selects how the mean shear (c_best) and its uncertainty
// (sigma) are estimated per spatial bin.
//   PDF_SIGMA          - chi2 sign test + quadratic fitting (c and sigma from PDF)
//   PDF_JACK - PDF chi2 sign test for c_best, jackknife for sigma
//   SWSE_JACK           - SWSE ratio estimator for c, jackknife for sigma
enum class StaticMode {
    PDF_SIGMA,  // PDF chi-square estimate with analytic sigma.
    PDF_JACK,   // PDF chi-square estimate with jackknife sigma.
    SWSE_JACK   // Ratio estimate with jackknife sigma.
};
inline constexpr StaticMode FD_STATIC_MODE = StaticMode::PDF_SIGMA;  // Active FD estimator mode.

// Derived helpers for compile-time branching
inline constexpr bool FD_USE_PDF_STATIS =   // Mode 1 or 2: statis uses PDF chi2 sign test
    (FD_STATIC_MODE == StaticMode::PDF_SIGMA ||
     FD_STATIC_MODE == StaticMode::PDF_JACK);
inline constexpr bool FD_USE_JACKKNIFE =    // Mode 2 or 3: plotComparison does jackknife for sigma
    (FD_STATIC_MODE == StaticMode::PDF_JACK ||
     FD_STATIC_MODE == StaticMode::SWSE_JACK);
inline constexpr bool FD_USE_SWSE_DATA =    // Mode 3: star cut uses SWSE data model
    (FD_STATIC_MODE == StaticMode::SWSE_JACK);

// Star-bar mode: true = per-exposure star bar (NtoN style),
//                false = single global star bar (Nto1 style)
inline constexpr bool FD_PER_EXPOSURE_STAR_BAR = false;  // Fit one star bar per exposure when true.

// ==================== Dimensions ====================
inline constexpr int nmax_per_core = 20000000;  // Maximum sources reserved per MPI rank.
inline constexpr int fd_num = 21;          // spatial bins by field distortion
inline constexpr int PDF_BINS = 4;         // equal-probability inner bins
inline constexpr float gf_lim = 0.0015;      // spatial bin range ±gf_lim
inline constexpr int NMAX = 200;           // fine grid sampling points
inline constexpr int MAX_DUP = 5;          // max duplicate measurements

// ==================== Jackknife / K-means ====================
inline constexpr int N_jack = 50;          // jackknife regions
inline constexpr int nmax_total = 1000000; // max total sources for k-means
inline constexpr int Km_iter = 100;        // k-means iterations

// ==================== Quality-cut thresholds ====================
inline constexpr float snrfcut = 4.0;  // Minimum Fourier signal-to-noise ratio.
inline constexpr float snrlow = 0.0;  // Optional lower source-SNR bound; zero disables it.
inline constexpr float snrhigh = 0.0;  // Optional upper source-SNR bound; zero disables it.
inline constexpr float starcut = 20.0;  // Point-source size cut.
inline constexpr float chi2_thresh = 0.01;  // Maximum exposure chi-square.
inline constexpr float flagcut = 0.0;  // Maximum accepted source quality flag.
inline constexpr float imaxcut = 64.0;  // Maximum source peak x coordinate.
inline constexpr float jmaxcut = 64.0;  // Maximum source peak y coordinate.
inline constexpr float zplow = 0.0;  // Minimum photometric redshift.
inline constexpr float zphigh = 3.0;  // Maximum photometric redshift.
inline constexpr float r_half_thresh = 0.0;  // Optional half-light-radius threshold.
inline constexpr float star_bar_mltp = 3.0;  // Stellar-locus sigma multiplier.
inline constexpr float psf_chi2_mltp = 3.0;  // PSF chi-square sigma multiplier.

// ==================== External-catalog cuts ====================
inline constexpr float ft_cut = -1.0;  // FLAGS_FT cut; negative disables it.
inline constexpr float fg_cut = -10.0;  // FLAGS_FG cut.
inline constexpr float gold_cut = -10.0;  // FLAGS_GOLD cut.
inline constexpr float ext_cut = -4.0;  // EXT_MASH cut.

// ==================== Star-cut histogram parameters ====================
inline constexpr int n_size_bins = 100;  // Stellar-size histogram bins.
inline constexpr int n_mag_bins = 20;  // Magnitude histogram bins.
inline constexpr float size_min = -2.0;  // Minimum histogram size coordinate.
inline constexpr float size_max = 2.0;  // Maximum histogram size coordinate.
inline constexpr float mag_min_val = 10.0;  // Minimum histogram magnitude.
inline constexpr float mag_max_val = 30.0;  // Maximum histogram magnitude.
inline constexpr int min_bin_count = 100;  // Minimum samples in a usable bin.
inline constexpr float peak_match_tol = 0.05;  // Stellar-peak matching tolerance.
inline constexpr float min_concentration = 0.6;  // Minimum stellar-locus concentration.
inline constexpr float star_phy_min = -0.5;  // Minimum physical stellar size.
inline constexpr float star_phy_max = 0.2;  // Maximum physical stellar size.

// Per-exposure star-cut additional parameters
inline constexpr float stage1_snr = 40.0;  // SNR for histogram accumulation
inline constexpr float stage2_snr = 0.0;  // Second-pass SNR threshold; zero disables it.
inline constexpr float init_win_active = 0.1;  // Active exposure initial size window.
inline constexpr float init_win_fallback = 0.15;  // Fallback exposure initial size window.
inline constexpr float default_s_init = 0.5;  // Default initial stellar-size center.
inline constexpr float clip_nsigma = 3.0;  // Iterative clipping sigma.
inline constexpr float min_clip_limit = 0.015;  // Minimum clipping half-width.
inline constexpr float default_s_std = 0.05;  // Default stellar-size scatter.
inline constexpr float fallback_scut_default = 0.6;  // Fallback stellar-size cut.

// ==================== Catalog column indices (0-based, DES format) ====================
// External-catalog columns (1-based in Fortran → 0-based here)
inline constexpr int col_flags_ft = 0;  // FLAGS_FT column index.
inline constexpr int col_flags_fg = 1;  // FLAGS_FG column index.
inline constexpr int col_flags_gold = 2;  // FLAGS_GOLD column index.
inline constexpr int col_ext_mash = 3;  // EXT_MASH column index.
inline constexpr int col_cra = 4;  // Right-ascension column index.
inline constexpr int col_cdec = 5;  // Declination column index.
inline constexpr int col_mag_g = 6;  // g-band magnitude column index.
inline constexpr int col_mag_r = 8;  // r-band magnitude column index.
inline constexpr int col_mag_i = 10;  // i-band magnitude column index.
inline constexpr int col_mag_z = 12;  // z-band magnitude column index.
inline constexpr int col_mag_y = 14;  // y-band magnitude column index.
inline constexpr int col_zp = 16;  // Photometric-redshift column index.
inline constexpr int external_num_cols =
    ExtCatConfig::EXTCAT_TOTAL_COLUMNS * LensingConfig::ext_cat;  // Optional external prefix width.
inline constexpr int col_expo = external_num_cols;  // Original 1-based exposure number.
inline constexpr int col_ccd = external_num_cols + 1;  // CCD number column index.
inline constexpr int source_col_offset = external_num_cols + 2;  // Fixed EXPO_NUM/CCD_NUM prefix end.

// Per-source columns (0-based absolute, derived from LensingConfig indices)
inline constexpr int col_polychi2 = source_col_offset + LensingConfig::iid;  // PSF fit chi-square column.
inline constexpr int col_pixx = source_col_offset + LensingConfig::ipixx;  // Source x column.
inline constexpr int col_pixy = source_col_offset + LensingConfig::ipixy;  // Source y column.
inline constexpr int col_sig = source_col_offset + LensingConfig::isig;  // Noise sigma column.
inline constexpr int col_star = source_col_offset + LensingConfig::istar;  // PSF-star count column.
inline constexpr int col_peak = source_col_offset + LensingConfig::ipeak;  // Historical peak column.
inline constexpr int col_imax = source_col_offset + LensingConfig::i_imax;  // Peak x column.
inline constexpr int col_jmax = source_col_offset + LensingConfig::i_jmax;  // Peak y column.
inline constexpr int col_h_flux = source_col_offset + LensingConfig::ih_flux;  // Half-light flux column.
inline constexpr int col_h_area = source_col_offset + LensingConfig::ih_area;  // Source area column.
inline constexpr int col_flag = source_col_offset + LensingConfig::iflag;  // Quality flag column.
inline constexpr int col_PSF = source_col_offset + LensingConfig::iPSF;  // Local PSF size column.
inline constexpr int col_SNR_F = source_col_offset + LensingConfig::iSNR_F;  // Fourier SNR column.
inline constexpr int col_ra = source_col_offset + LensingConfig::ira;  // Source RA column.
inline constexpr int col_dec = source_col_offset + LensingConfig::idec;  // Source Dec column.
inline constexpr int col_gf1 = source_col_offset + LensingConfig::igf1;  // Field distortion g1 column.
inline constexpr int col_gf2 = source_col_offset + LensingConfig::igf2;  // Field distortion g2 column.
inline constexpr int col_g1 = source_col_offset + LensingConfig::ig1;  // Fourier_Quad g1 column.
inline constexpr int col_g2 = source_col_offset + LensingConfig::ig2;  // Fourier_Quad g2 column.
inline constexpr int col_de = source_col_offset + LensingConfig::ide;  // Shear response column.
inline constexpr int col_h1 = source_col_offset + LensingConfig::ih1;  // Higher-order h1 column.
inline constexpr int col_h2 = source_col_offset + LensingConfig::ih2;  // Higher-order h2 column.
inline constexpr int col_cos2 = source_col_offset + LensingConfig::icos2;  // Spin-2 cosine column.
inline constexpr int col_sin2 = source_col_offset + LensingConfig::isin2;  // Spin-2 sine column.
inline constexpr int col_parity = source_col_offset + LensingConfig::iparity;  // WCS parity column.
inline constexpr int col_chi2 = source_col_offset + LensingConfig::ichi2;  // Exposure chi-square column.

// Total number of columns in the catalog
inline constexpr int ICHI2 = source_col_offset + LensingConfig::npara;  // Total catalog row width.

static_assert(col_expo + 1 == col_ccd,
              "EXPO_NUM must immediately precede CCD_NUM");
static_assert(col_ccd + 1 == col_polychi2,
              "Stage-7 fields must immediately follow CCD_NUM");
static_assert(col_chi2 + 1 == ICHI2,
              "Chi2 must be the final catalog column");

// ==================== Bad CCD list (DES) ====================
inline constexpr int bad_ccds[] = {2, 31, 53, 61};  // DES CCD numbers excluded from analysis.
inline constexpr int n_bad_ccds = 4;  // Number of excluded CCDs.

// ==================== Chip-edge masking (DES) ====================
inline constexpr int chip_xmin = 50;  // Minimum accepted chip x coordinate.
inline constexpr int chip_xmax = 1990;  // Maximum accepted chip x coordinate.
inline constexpr int chip_ymin = 100;  // Minimum accepted chip y coordinate.
inline constexpr int chip_ymax = 3990;  // Maximum accepted chip y coordinate.

}  // namespace FDConfig

#endif  // FD_CONFIG_HPP
