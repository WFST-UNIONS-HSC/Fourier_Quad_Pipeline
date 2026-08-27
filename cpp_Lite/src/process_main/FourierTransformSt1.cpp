#include "process_main/FourierTransformSt1.hpp"
#include "process_main/ProcessMainState.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/ImageProcessing.hpp"
#include "process_main/MPIFailure.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

namespace FourierTransformSt1 {

    // ==========================================
    // Function: Run Stage 4 for every chip in one Lite exposure
    // Method: Resolve the exposure list and delegate each chip to the norm-gated transformer.
    // ==========================================
    void procFourierTSt1(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
        std::vector<std::string> image_files;
        std::string dir_output;
        UniversalUtils::getImageList(expo_file_path, image_files, dir_output);

        for (const auto& image_file : image_files) {
            chipProcessFourierTSt1(image_file, dir_output);
        }
    }

    // ==========================================
    // Function: Transform one chip's star-candidate stamps to Fourier power
    // Method: Prepare the configured noise product, build one shared corrected-power path,
    //         regularize it, and size all buffers from the live candidate count.
    // ==========================================
    void chipProcessFourierTSt1(const std::string& imageFile, const std::string& dirOutput) {
        const Universalblock::NormStatus norm_status =
            Universalblock::checkNorm(imageFile, dirOutput);
        if (norm_status == Universalblock::NormStatus::Invalid) {
            return;
        }
        if (norm_status != Universalblock::NormStatus::Valid) {
            MPIFailure::abortWorld(
                "check Stage 1 norm before star FFT",
                Universalblock::normErrorDetail(
                    norm_status, imageFile, dirOutput));
        }

        std::string raw_prefix = UniversalUtils::getPrefix(imageFile);
        // PREFIX inlined: per-type stamps/ subdirs (reorganized layout)

        int nsource = 0;
        std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/dat_StarCanInfo", raw_prefix, "_star_can_info.dat");

        std::ifstream fin(filename);
        if (!fin.is_open()) {
            MPIFailure::abortWorld("read star-candidate info", filename);
        }

        std::string header;
        std::getline(fin, header); // skip header line

        double val1, val2, val3, val4;
        while (fin >> val1 >> val2 >> val3 >> val4) {
            nsource++;
        }
        fin.close();

        if (nsource > 0) {
            int ns = LensingConfig::ns;
            int len_s = LensingConfig::len_s;
            int nn1 = ns * len_s;
            int nn2 = ns * (nsource / len_s + 1);

            std::vector<float> source_coll;
            std::vector<float> noise_coll;

            std::string filename_star_can = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCan", raw_prefix, "_star_can.fits");
            if (!FitsIO::readStamps(nsource, 1, nsource, ns, ns, source_coll, nn1, nn2, filename_star_can)) {
                MPIFailure::abortWorld(
                    "read star-candidate stamps", filename_star_can);
            }

            std::string filename_star_can_noise = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanN", raw_prefix, "_star_can_noise.fits");
            if (!FitsIO::readStamps(nsource, 1, nsource, ns, ns, noise_coll, nn1, nn2, filename_star_can_noise)) {
                MPIFailure::abortWorld(
                    "read star-candidate noise stamps",
                    filename_star_can_noise);
            }

            std::vector<float> power_coll(
                static_cast<size_t>(nsource) * ns * ns, 0.0f);

            for (int i = 0; i < nsource; ++i) {
                std::vector<float> source(ns * ns);
                std::vector<float> noise_product(ns * ns);
                std::copy(source_coll.begin() + i * ns * ns,
                          source_coll.begin() + (i + 1) * ns * ns,
                          source.begin());
                std::copy(noise_coll.begin() + i * ns * ns,
                          noise_coll.begin() + (i + 1) * ns * ns,
                          noise_product.begin());

                std::vector<float> source_p(ns * ns);
                std::vector<float> noise_p;
                double source_pc = 0.0;

                if (!ImageProcessing::prepareNoisePower(
                        ns, noise_product, LensingConfig::NstampType, noise_p)) {
                    MPIFailure::abortWorld(
                        "prepare star-candidate noise power", filename_star_can_noise);
                }
                if (!ImageProcessing::buildCorrectedPower(
                        ns, ns, source, noise_p, LensingConfig::star_smooth,
                        source_p, source_pc)) {
                    MPIFailure::abortWorld(
                        "build corrected star-candidate power", filename_star_can);
                }
                ImageProcessing::regularizePower(ns, ns, source_p, LensingConfig::star_smooth);

                std::copy(source_p.begin(), source_p.end(), power_coll.begin() + i * ns * ns);
            }

            std::string filename_star_can_power = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", raw_prefix, "_star_can_power.fits");
            FitsIO::writeStamps(
                nsource, 1, nsource, ns, ns, power_coll, nn1, nn2,
                filename_star_can_power);
        }
    }
}
