#include "LensingConfig.hpp"
#include "MPIScheduler.hpp"
#include "PreProcess.hpp"
#include "Astrometry.hpp"
#include "SourceExtractor.hpp"
#include "FourierTransformSt1.hpp"
#include "FourierTransformSt2.hpp"
#include "PSFModel.hpp"
#include "PSFRecons.hpp"
#include "ShearMeasurement.hpp"
#include "ExposureInfo.hpp"
#include "CatalogCombiner.hpp"
#include "UniversalUtils.hpp"
#include "NumericalRecipes.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

std::vector<std::string> EXPO_FILE;
int N_EXPO = 0;

void initialize(const std::string& expoList) {
    EXPO_FILE.clear();
    std::ifstream fin(expoList);
    if (!fin.is_open()) {
        std::cerr << "EXPO_LIST reading error: " << expoList << std::endl;
        std::exit(1);
    }
    std::string expo_name;
    int nchip;
    while (fin >> expo_name >> nchip) {
        // EXPO_LIST Format : "DIR_OUTPUT/stamps/prefix_expo.list"     NChip
        if (expo_name.size() >= 2 && expo_name.front() == '"' && expo_name.back() == '"') {
            expo_name = expo_name.substr(1, expo_name.size() - 2);
        }
        EXPO_FILE.push_back(expo_name);
    }
    fin.close();
    N_EXPO = static_cast<int>(EXPO_FILE.size());
    std::cout << "Total number of EXPOSURE: " << N_EXPO << std::endl;
}

// ==========================================
// Function: Run the MPI Fourier_Quad pipeline
// Method: Initialize MPI, seed each rank's F77-equivalent RNG once, then execute enabled stages.
// ==========================================
int main(int argc, char* argv[]) {
    // 1. MPI Initialization
    MPIScheduler::init(argc, argv);
    int my_id = MPIScheduler::my_id;

    // Match F77: read each process clock, mix it with rank, and seed ran1 exactly once.
    const unsigned int rng_seed = NumericalRecipes::initializeRan1Seed(
        my_id, MPIScheduler::num_procs);
    std::cout << "RNG_SEED rank seed: " << my_id << " " << rng_seed << std::endl;
    MPIScheduler::barrier();

    if (argc < 2) {
        if (my_id == 0) {
            std::cerr << "Usage: " << argv[0] << " <EXPO_LIST>" << std::endl;
        }
        MPIScheduler::finalize();
        return 1;
    }
    std::string EXPO_LIST = argv[1];

    // 2. Initialize exposure list on Rank 0 and broadcast to all ranks
    if (my_id == 0) {
        initialize(EXPO_LIST);
    }
    
    // Broadcast N_EXPO
    MPI_Bcast(&N_EXPO, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize EXPO_FILE on non-zero ranks
    if (my_id != 0) {
        EXPO_FILE.resize(N_EXPO);
    }

    // Broadcast string lengths and contents
    for (int i = 0; i < N_EXPO; ++i) {
        int str_len = 0;
        if (my_id == 0) {
            str_len = static_cast<int>(EXPO_FILE[i].length());
        }
        MPI_Bcast(&str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (my_id != 0) {
            EXPO_FILE[i].resize(str_len);
        }
        MPI_Bcast(const_cast<char*>(EXPO_FILE[i].data()), str_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    MPIScheduler::barrier();

    // ==========================================
    // Function: Validate stage dependency before execution
    // Method: Stage 9 consumes Stage 8 exposure chi2, so reject PROCESS_stage with 23 but without 19.
    // ==========================================
    if (LensingConfig::PROCESS_stage % 23 == 0 && LensingConfig::PROCESS_stage % 19 != 0) {
        if (my_id == 0) {
            std::cerr << "Error: Stage 9 requires Stage 8. PROCESS_stage enables CatalogCombiner without ExposureInfo." << std::endl;
        }
        MPIScheduler::finalize();
        return 1;
    }

    // 3. Stage 1: Pre-process (mod(PROCESS_stage, 2) == 0)
    if (LensingConfig::PROCESS_stage % 2 == 0) {
        MPIScheduler::distribute(N_EXPO, PreProcess::preProcess, "Pre-process...");
    }
    MPIScheduler::barrier();

    // 4. Stage 2: Astrometry (mod(PROCESS_stage, 3) == 0)
    if (LensingConfig::PROCESS_stage % 3 == 0) {
        MPIScheduler::distribute(N_EXPO, Astrometry::procAstrometry, "Astrometry...");
    }
    MPIScheduler::barrier();

    // 5. Stage 3: Sources (mod(PROCESS_stage, 5) == 0)
    if (LensingConfig::PROCESS_stage % 5 == 0) {
        MPIScheduler::distribute(N_EXPO, SourceExtractor::procSource, "Sources ...");
    }
    MPIScheduler::barrier();

    // 6. Stage 4: FFT st1 (mod(PROCESS_stage, 7) == 0)
    if (LensingConfig::PROCESS_stage % 7 == 0) {
        MPIScheduler::distribute(N_EXPO, FourierTransformSt1::procFourierTSt1, "FFT st1...");
    }
    MPIScheduler::barrier();

    // 7. Stage 5: PSF (mod(PROCESS_stage, 11) == 0)
    if (LensingConfig::PROCESS_stage % 11 == 0) {
        MPIScheduler::distribute(N_EXPO, PSFModel::procPSF, "PSF ...");
        if (LensingConfig::PSF_Ms == 1) {
            PSFRecons::chipPSFRecons(N_EXPO);
        }
    }
    MPIScheduler::barrier();

    // 8. Stage 6: FFT st2 (mod(PROCESS_stage, 13) == 0)
    if (LensingConfig::PROCESS_stage % 13 == 0) {
        MPIScheduler::distribute(N_EXPO, FourierTransformSt2::procFourierTSt2, "FFT st2...");
    }
    MPIScheduler::barrier();

    // 9. Stage 7: Shear (mod(PROCESS_stage, 17) == 0)
    if (LensingConfig::PROCESS_stage % 17 == 0) {
        MPIScheduler::distribute(N_EXPO, ShearMeasurement::procShear, "Shear ...");
    }
    MPIScheduler::barrier();

    if (LensingConfig::PSF_Ms == 1) {
        PSFModel::freePSFMemory();
    }

    // 10. Stage 8: Info (mod(PROCESS_stage, 19) == 0)
    ExposureInfo::expo_para.assign(static_cast<size_t>(LensingConfig::NMAX_EXPO) * 6, 0.0f);
    if (LensingConfig::PROCESS_stage % 19 == 0) {
        MPIScheduler::distribute(N_EXPO, ExposureInfo::procInfo, "Info ...");
    }
    MPIScheduler::barrier();

    // Gather exposure parameters across all processes
    std::vector<float> expo_para_t(static_cast<size_t>(LensingConfig::NMAX_EXPO) * 6, 0.0f);
    MPI_Allreduce(ExposureInfo::expo_para.data(), expo_para_t.data(),
                  LensingConfig::NMAX_EXPO * 6, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    ExposureInfo::expo_para = expo_para_t;

    if (my_id == 0) {
        std::string root_dir = UniversalUtils::getDir(EXPO_LIST, 1);
        std::string filename = root_dir + "/expo_info.dat";
        std::ofstream fout(filename);
        if (fout.is_open()) {
            fout << std::setprecision(10);
            fout << "N-valid-chip PSF-FWHM(arcsec) chi_d-stars nstar-per-chip cRVAL1 cRVAL2 expo_name\n";
            for (int i = 0; i < N_EXPO; ++i) {
                for (int j = 0; j < 6; ++j) {
                    fout << ExposureInfo::expo_para[i * 6 + j] << " ";
                }
                fout << EXPO_FILE[i] << "\n";
            }
            fout.close();
        } else {
            std::cerr << "Error: cannot write to expo_info.dat file: " << filename << std::endl;
        }
    }
    MPIScheduler::barrier();

    // 11. Stage 9: Combine (mod(PROCESS_stage, 23) == 0)
    if (LensingConfig::PROCESS_stage % 23 == 0) {
        MPIScheduler::distribute(N_EXPO, CatalogCombiner::procComb, "combine ...");
    }
    MPIScheduler::barrier();

    // 12. MPI Finalization
    MPIScheduler::barrier();
    MPIScheduler::finalize();

    return 0;
}
