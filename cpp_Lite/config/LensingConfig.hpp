#ifndef LENSING_CONFIG_HPP
#define LENSING_CONFIG_HPP

#include "Initialize.hpp"
#include "pathconfig.hpp"

#include <cmath>
// ==========================================
// cpp_lite: the following build-time branch selectors of the full pipeline are FROZEN and
// their unused branches have been deleted from the sources.  They no longer exist as
// constants; the code implements exactly one behaviour each:
//   ASTROMETRY_trivial = 0  -> Gaia-based astrometry only
//   include_FLAT       = 0  -> no super-flat multiplication
//   include_Mask       = 2  -> per-chip DQ mask from dirOutput/dqmask
//   PreprocsType       = 1  -> historical F77 preprocessing
//   include_BGsub      = 0  -> no Stage-3 background re-subtraction
//   ext_cat            = 1  -> configured external source catalogue
//   ext_PSF            = 0  -> PSF measured from the stars of the frame
//   deblending         = 1  -> de-blending always applied
//   NstampType         = 1  -> deterministic historical F77 blank noise
//   PsfGroupingType    = 1  -> historical F77 largest-component grouping
//   PSF_type           = 1  -> local polynomial PSF fit
//   PSF_Ms             = 0  -> no multi-scale / PCA PSF reconstruction
// Still selectable: PROCESS_stage, CCD_split, gal_smooth, and star_smooth.
// ==========================================
namespace LensingConfig {
    // Camera geometry
    inline constexpr int N_CCD = Initialize::N_CCD;  // Number of CCD chips per exposure.
    // CCD configuration
    inline constexpr double pixel_size = Initialize::pixel_size;  // Pixel scale in arcseconds.
    inline constexpr double saturation_thresh = Initialize::saturation_thresh;  // Saturated pixel threshold.
    // Image/CCD size parameters
    inline constexpr int chipnx = Initialize::chipnx;  // Science CCD width used for PSF coordinates.
    inline constexpr int chipny = Initialize::chipny;  // Science CCD height used for PSF coordinates.
    // Stage control parameters
    inline constexpr int PROCESS_stage = Initialize::PROCESS_stage;  // Prime-product stage selector.

    // Split parameters
    inline constexpr int CCD_split = Initialize::CCD_split;  // Split each CCD into amplifier regions.
    constexpr int nct = 12;  // Number of background rectangles.
    constexpr int ncx = 3;  // Number of background rectangles along x.

    // PSF selection and configuration
    constexpr int psf_order = 8;  // PSF polynomial order selector.
    constexpr int npo = 64;  // Exposure PSF sample count.
    constexpr int npox = 8;  // Exposure PSF samples along x.
    constexpr int nstar_min = npo * 3 / 2;  // Minimum stars for exposure PSF fitting.
    constexpr int npl = 10;  // Local PSF polynomial coefficient count minus one.
    constexpr int nplx = 2;  // Local PSF polynomial degree along x.
    constexpr int nstar_min_local = 16;  // Minimum retained stars for a local fit.

    // ==========================================
    // Configuration: Frozen Stage-5 initial local-PSF fit
    // Method: Retain the leverage denominator guard used by the one live fit.
    // ==========================================
    constexpr double psf_loo_min_denom = 1.0e-6;  // Minimum leave-one-out denominator.
    // ==========================================

    // Stamp dimensions
    constexpr int ns = 64;  // Science stamp and Fourier-grid side length.
    constexpr int nsns = ns * ns;  // Pixels in one science stamp.
    constexpr int chip_margin = 8;  // Extra chip-edge extraction margin.
    constexpr int ns_2 = ns / 2;  // Half science-stamp side length.
    constexpr int nl_2 = ns_2 + chip_margin;  // Half expanded extraction side.
    constexpr int nl = nl_2 * 2;  // Full expanded extraction side.
    constexpr int flag_thresh = 3;  // Maximum accepted source extraction flag.
    constexpr int chip_edge_margin = chip_margin;  // Alias used by chip-edge checks.

    constexpr double dz_thresh = 0.1;  // Redshift tolerance for catalog matching.

    // Catalog sizes and limits
    constexpr int len_g = 40;  // Galaxy metadata row capacity.
    constexpr int len_s = 15;  // Star metadata row capacity.
    // Maximum number of flux-ranked image detections passed to astrometric pattern matching.
    // This is a scientific selection limit, not a catalog-storage capacity limit.
    constexpr int n_user_max = 200;  // Bright detections used for astrometric matching.
    constexpr int ngal_max = 2000;  // Initial galaxy-vector reservation hint.
    constexpr int nstar_max = 1000;  // Initial star-vector reservation hint.
    constexpr int npara = 25;  // Per-source Stage-7 catalog field count.
    constexpr int len_sam = 50;  // PSF sample metadata row length.

    constexpr int npd = 33;  // PU astrometric distortion coefficient count.
    // Target side length for the balanced background blocks.
    constexpr int blocksize = 200;  // Target background block side length.

    // Thresholds
    constexpr double source_thresh = 2.0;  // Source-detection SNR threshold.
    constexpr double core_thresh = 4.0;  // Source-core detection threshold.

    // ==========================================
    // Configuration: Retained local noise-plane fitting utility
    // Method: Define the square shell and minimum valid fraction exercised by
    //         the standalone NoisePlaneFit regression.
    // ==========================================
    constexpr int noise_region_size = 192;  // Outer local-noise square side length.
    constexpr int noise_inner_size = 96;  // Central exclusion square side length.
    constexpr double noise_plane_min_valid_fraction = 0.30;  // Minimum plane-fit shell fraction.

    constexpr int area_max = ns * ns;  // Maximum connected source area.
    constexpr int area_thresh = 6;  // Minimum connected source area.

    constexpr int gal_smooth = 0;  // Galaxy-stamp smoothing type.
    constexpr int star_smooth = 2;  // Star-stamp smoothing type.
    constexpr double SNR_PSF = 100.0;  // Minimum PSF-star signal-to-noise ratio.

    // Catalogue column indices (shifted to 0-based for C++)
    constexpr int iid = 1 - 1;  // PSF polynomial chi-square field index.
    constexpr int ipixx = 2 - 1;  // Source-center x field index.
    constexpr int ipixy = 3 - 1;  // Source-center y field index.
    constexpr int isig = 4 - 1;  // Local noise sigma field index.
    constexpr int istar = 5 - 1;  // Available PSF-star count field index.
    constexpr int ipeak = 5 - 1;  // Historical peak alias field index.
    constexpr int i_imax = 6 - 1;  // Peak x field index.
    constexpr int i_jmax = 7 - 1;  // Peak y field index.
    constexpr int ih_flux = 8 - 1;  // Half-light flux field index.
    constexpr int ih_area = 9 - 1;  // Source area field index.
    constexpr int iflag = 10 - 1;  // Quality flag field index.
    constexpr int iPSF = 11 - 1;  // Local PSF size field index.
    constexpr int iSNR_F = 12 - 1;  // Fourier SNR field index.
    constexpr int ira = 13 - 1;  // Right-ascension field index.
    constexpr int idec = 14 - 1;  // Declination field index.
    constexpr int igf1 = 15 - 1;  // Field-distortion g1 index.
    constexpr int igf2 = 16 - 1;  // Field-distortion g2 index.
    constexpr int ig1 = 17 - 1;  // Fourier_Quad g1 estimator index.
    constexpr int ig2 = 18 - 1;  // Fourier_Quad g2 estimator index.
    constexpr int ide = 19 - 1;  // Shear response estimator index.
    constexpr int ih1 = 20 - 1;  // Higher-order h1 estimator index.
    constexpr int ih2 = 21 - 1;  // Higher-order h2 estimator index.
    constexpr int icos2 = 22 - 1;  // Spin-2 cosine field index.
    constexpr int isin2 = 23 - 1;  // Spin-2 sine field index.
    constexpr int iparity = 24 - 1;  // WCS parity field index.
    constexpr int ichi2 = 25 - 1;  // Exposure chi-square field index.

    // Band correction parameters
    constexpr double g1_c = 0.0;  // Additive field-distortion g1 correction.
    constexpr double g2_c = 0.0;  // Additive field-distortion g2 correction.
    constexpr double chi2_thresh = 0.1;  // Maximum exposure PSF chi-square.

    // Mathematical constants
    constexpr double pi = 3.14159265358979323846;  // Mathematical pi.
    constexpr double arc_convert = pi / 180.0;  // Degrees-to-radians conversion factor.

    // ==========================================
    static_assert(psf_loo_min_denom > 0.0 && psf_loo_min_denom < 1.0,
                  "PSF LOO denominator floor must lie in (0,1)");
    static_assert(noise_region_size > noise_inner_size,
                  "noise region must exceed the central exclusion");
    static_assert(noise_inner_size >= nl,
                  "noise inner exclusion must cover the full source extraction region");
    static_assert(noise_region_size % 2 == 0 && noise_inner_size % 2 == 0,
                  "noise region and exclusion sizes must be even");
    static_assert((noise_region_size - noise_inner_size) % 2 == 0,
                  "noise inner exclusion must be centered on the local noise region");
    static_assert(noise_plane_min_valid_fraction > 0.0
                      && noise_plane_min_valid_fraction <= 1.0,
                  "noise plane minimum valid fraction must lie in (0,1]");
}

#endif // LENSING_CONFIG_HPP
