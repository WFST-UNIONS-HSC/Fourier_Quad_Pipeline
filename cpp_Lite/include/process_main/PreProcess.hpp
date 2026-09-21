#ifndef PRE_PROCESS_HPP
#define PRE_PROCESS_HPP

#include <string>
#include <vector>
#include "process_main/FitsIO.hpp"

namespace PreProcess {

    // Stage 1 driver
    void preProcess(int iexpo);

    // Individual chip preprocessing
    void chipPreProcess(const std::string& imageFile, const std::string& dirOutput, int cid);

    // Helper functions for preprocessing
    // ==========================================
    // Function: Estimate and subtract one amplifier's background model
    // Method: Run the frozen historical F77 estimator and publish coefficients in
    //         the normalized-coordinate FITS metadata basis.
    // ==========================================
    void setBackground(int x_start, int x_end, int y_start, int y_end, int nx, int ny,
                       std::vector<float>& image, int blocksize, int nct, int ncx,
                       std::vector<double>& bg_coeffs,
                       int& ierror);

    // ==========================================
    // Function: Estimate and apply one amplifier's historical noise-sigma plane
    // Method: Run the frozen random-2000 F77 estimator with its positivity guard.
    // ==========================================
    void setSig(int x_start, int x_end, int y_start, int y_end, int nx, int ny, std::vector<float>& image,
                double& aa, double& bb, double& cc, int& ierror);

    void locateDefects(int nx, int ny, const std::vector<float>& array, std::vector<float>& normap,
                       std::vector<int>& weight, int area_max, int area_thresh, int& ierror);

    void mergeDefects(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                      int area_max, double source_thresh, int area_thresh, int& ierror);

    void detectArtificialStripes(int nx, int ny, std::vector<int>& weight,
                                 const std::vector<float>& diffx, const std::vector<float>& diffy,
                                 float sigx, float sigy, float medx, float medy);

    void maskSourceRegions(int nx, int ny, std::vector<int>& weight, const std::vector<float>& normap,
                           int area_max, double source_thresh, int area_thresh);

    void detectStripes(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                       int x_smooth, int y_smooth);

    void detectStellarHalo(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                           int npmax, double defect_halo_thresh);

    void detectDent(int nx, int ny, const std::vector<float>& normap, std::vector<int>& weight,
                    int npmax, double defect_halo_thresh);

}

#endif // PRE_PROCESS_HPP
