#include "FourierTransformSt2.hpp"
#include "OutputFile.hpp"
#include "OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"
#include "Universalblock.hpp"
#include "FitsIO.hpp"
#include "ImageProcessing.hpp"
#include "NoiseCovariance.hpp"
#include "MPIFailure.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>

extern std::vector<std::string> EXPO_FILE;

namespace FourierTransformSt2 {

// ==========================================
// Function: Transform one chip's source stamps into Fourier-space products
// Method: Apply the shared norm gate, size buffers from live rows, update diagnostics,
//         and publish text/FITS outputs through checked main-process writers.
// ==========================================
void chipProcessFourierTSt2(const std::string& imageFile, const std::string& dirOutput) {
    const Universalblock::NormStatus norm_status =
        Universalblock::checkNorm(imageFile, dirOutput);
    if (norm_status == Universalblock::NormStatus::Invalid) {
        return;
    }
    if (norm_status != Universalblock::NormStatus::Valid) {
        MPIFailure::abortWorld(
            "check Stage 1 norm before galaxy FFT",
            Universalblock::normErrorDetail(norm_status, imageFile, dirOutput));
    }

    int ns = LensingConfig::ns;

    std::string raw_prefix = UniversalUtils::getPrefix(imageFile);
    // PREFIX inlined: per-type stamps/ subdirs (reorganized layout)

    int nsource = 0;
    std::string info_filename = OutputLayout::chipPath(
        dirOutput, "stamps/dat_SrcInfo", raw_prefix, "_source_info.dat");
    std::vector<std::vector<float>> source_para;

    std::ifstream fin(info_filename);
    if (!fin.is_open()) {
        MPIFailure::abortWorld("read source info for galaxy FFT", info_filename);
    }

    std::string header;
    std::getline(fin, header); // skip header line

    std::string line;
    int iflag_col = LensingConfig::iflag + 1; // 10 columns
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<float> row(LensingConfig::npara, 0.0f);
        bool success = true;
        for (int i = 0; i < iflag_col; ++i) {
            if (!(ss >> row[i])) {
                success = false;
                break;
            }
        }
        if (success) {
            source_para.push_back(row);
        }
    }
    fin.close();

    nsource = source_para.size();
    if (nsource == 0) {
        // Write header and return
        MainIO::OutputFile fout(info_filename);
        fout << "ig xp yp sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F\n";
        fout.close();
        return;
    }

    int len_g = LensingConfig::len_g;
    int nn1 = ns * len_g;
    int nn2 = ns * (nsource / len_g + 1);

    std::vector<float> source_coll;
    std::string source_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Src", raw_prefix, "_source.fits");
    if (!FitsIO::readStamps(nsource, 1, nsource, ns, ns, source_coll, nn1, nn2, source_fits)) {
        MPIFailure::abortWorld("read source stamps for galaxy FFT", source_fits);
    }

    std::vector<float> noise_coll;
    std::string noise_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_Noise", raw_prefix, "_noise.fits");
    if (!FitsIO::readStamps(nsource, 1, nsource, ns, ns, noise_coll, nn1, nn2, noise_fits)) {
        MPIFailure::abortWorld("read noise stamps for galaxy FFT", noise_fits);
    }

    std::vector<float> power_coll(
        static_cast<size_t>(nsource) * ns * ns, 0.0f);

    for (int i = 0; i < nsource; ++i) {
        std::vector<float> source(ns * ns);
        std::copy(source_coll.begin() + i * ns * ns, source_coll.begin() + (i + 1) * ns * ns, source.begin());

        std::vector<float> source_p(ns * ns);
        std::vector<float> noise_p;
        double pc = 0.0;

        // SNR calculation uses star_smooth (2)
        ImageProcessing::getPower(ns, ns, source, source_p, 2, pc);

        int ns_2 = LensingConfig::ns_2;
        float cen_val = source_p[ns_2 * ns + ns_2];
        source_para[i][10] = std::sqrt(std::max(static_cast<float>(pc), cen_val));
        source_para[i][11] = source_para[i][10] / source_para[i][3] * ns;

        // Main power spectrum computation uses gal_smooth
        ImageProcessing::getPower(ns, ns, source, source_p, LensingConfig::gal_smooth, pc);
        if (!NoiseCovariance::copyStoredNoisePower(
                noise_coll, static_cast<std::size_t>(i) * ns * ns, ns, noise_p)) {
            MPIFailure::abortWorld("load stored galaxy noise power", noise_fits);
        }
        ImageProcessing::processPowers(ns, source_p, noise_p);

        std::copy(source_p.begin(), source_p.end(), power_coll.begin() + i * ns * ns);
    }

    MainIO::OutputFile fout(info_filename);
    if (!fout.is_open()) {
        std::cerr << "Error opening " << info_filename << " for output" << std::endl;
        return;
    }
    fout << "ig xp yp sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F\n";
    for (int i = 0; i < nsource; ++i) {
        for (int j = 0; j <= LensingConfig::iSNR_F; ++j) {
            fout << source_para[i][j] << (j == LensingConfig::iSNR_F ? "" : " ");
        }
        fout << "\n";
    }
    fout.close();

    std::string power_fits = OutputLayout::chipPath(
        dirOutput, "stamps/fits_SrcP", raw_prefix, "_source_p.fits");
    FitsIO::writeStamps(nsource, 1, nsource, ns, ns, power_coll, nn1, nn2, power_fits);
}

// ==========================================
// Function: Run Stage 6 for every chip in one exposure
// Method: Resolve the exposure list and delegate each chip to the norm-gated transformer.
// ==========================================
void procFourierTSt2(int iexpo) {
    if (iexpo <= 0 || iexpo > static_cast<int>(EXPO_FILE.size())) {
        std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
        return;
    }
    std::string expo_file_path = EXPO_FILE[iexpo - 1];
    std::vector<std::string> image_files;
    std::string dir_output;
    UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

    for (const auto& image_file : image_files) {
        chipProcessFourierTSt2(image_file, dir_output);
    }
}

} // namespace FourierTransformSt2
