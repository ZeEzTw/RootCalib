/**
 * @class Histogram
 * @brief Provides a comprehensive suite of histogram management operations,
 *        including calibration, peak detection, and data export.
 *
 *        This class offers the following core functionalities:
 *
 *        - **Peak Detection**: Detects peaks within specified ranges using an in-house
 *          `Peak` class, all peaks are fitted by Gaussian curves.
 *
 *        - **Spectrum Calibration**: Calibrates the spectrum across all defined degrees,
 *          delivering an adjusted spectrum based on user-specified calibration sources.
 *
 *
 *
 *        - **Data Export**: Outputs results in JSON format and ROOT files, making it easier
 *          to share and analyze calibrated data and peak characteristics.
 *
 *        All functions in this class are custom-built, including peak detection,
 *        calibration degree determination, and mathematical problem-solving utilities,
 *        utilizing only the standard capabilities of C++.
 *
 *        - **Extra Features**:
 *             -Peaks that are not used in calibration, but are detected are marked with a green
 *             gaussian in the histogram instead of a red one, better for debugging.
 *
 *        Additional Features: avable in previous versions of the code (on github)
 *
 *        - **Degree-Free Calibration**: Performs calibration without predefined thresholds
 *          to determine the most suitable polynomial degree.
 *
 *        - **ROOT-Based Calibration**: An alternative calibration method leveraging ROOT
 *          functions for users requiring integration with ROOT analysis tools.
 *
 *        - **High-Precision Polynomial Calibration**: Supports precise polynomial
 *          calibration up to degree 1, enhancing accuracy for users with detailed calibration needs.
 */

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Peak.h"
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>

class Histogram
{
private:
    // Core histogram data
    TH1D *mainHist;
    TH1D *tempHist;
    TH1D *calibratedHist;
    std::vector<Peak> peaks;
    std::vector<double> coefficients;

    // Configuration parameters
    int xMin, xMax;
    int maxFWHM;
    float minAmplitude;
    float maxAmplitude;
    int numberOfPeaks;
    float polynomialFitThreshold;

    // Calibration parameters
    int calibrationDegree;
    float m, b;
    unsigned int peakMatchCount;
    int peakCount;
    int polynomialDegree;

    // Metadata
    std::string TH2histogram_name;
    std::string sourceName;
    std::string serial;
    int detType;

    // Area calculations
    float totalArea;
    float totalAreaError;

    // Private methods for peak detection and fitting
    void eliminatePeak(const Peak &peak);
    TF1 *createGaussianFit(int maxBin);
    int findMaxBin();
    int detectAndFitPeaks();
    bool isValidPeak(const Peak &peak) const;
    bool checkConditions(const Peak &peak) const;
    bool extractDataForFit(std::vector<double> &xValues, std::vector<double> &yValues);
    bool areCoefficientsValid(TF1 *fitFunction, int degree, double threshold);
    void setCalibratedBinContent(int original_bin);
    std::string getMainHistName() const;
    void findStartOfPeak(Peak &peak);

    // Private methods for calibration
    double evaluateCalibrationPolynomial(double x) const;
    void initializeCalibratedHist();
    void interpolateBins(int start_bin, int end_bin, double start_value, double end_value,
                         double start_position, double end_position);
    bool checkPredictedEnergies(double predictedEnergy, const double knownEnergies[],
                                int size, float errorAdmitted, double &valueAssociatedWith) const;

public:
    // Constructors and destructor
    Histogram();
    Histogram(int xMin, int xMax, int maxFWHM, float minAmplitude, float maxAmplitude,
              const std::string &serial, int detType, float polynomialFitThreshold,
              int numberOfPeaks, TH1D *mainHist, const std::string &TH2histogram_name,
              const std::string &sourceName);
    ~Histogram();
    Histogram(const Histogram &histogram);
    Histogram &operator=(const Histogram &histogram);

    // Core functionality
    void findPeaks();
    void calibratePeaks(const double knownEnergies[], int size);
    void calibratePeaksByDegree();
    void applyXCalibration();
    void changePeak(int peakNumber, double newPosition);

    // Output methods
    void outputPeaksDataJson(std::ofstream &file);
    void printHistogramWithPeaksRoot(TFile *outputFile);
    void printCalibratedHistogramRoot(TFile *outputFile) const;

    // Getters and setters
    TH1D *getCalibratedHist() const { return calibratedHist; }
    TH1D *getMainHist() const { return mainHist; }
    unsigned int getPeakMatchCount() const { return peakMatchCount; }
    float getPT();
    float getPTError();
    void setTotalArea();
    void setTotalAreaError();
    const char *returnNameOfHistogram() const;
};

#endif // HISTOGRAM_H
