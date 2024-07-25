#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Peak.h"
#include <TH1D.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h> // Include corect pentru TLatex
class Histogram {
private:
    int xMin, xMax;
    int maxFWHM;
    int numberOfPeaks;
    TH1D *mainHist;
    TH1D *tempHist;
    TH1D *calibratedHist;
    float m, b;
    std::vector<Peak> peaks;
    int peakCount;

    // Funcții private
    void eliminatePeak(const Peak& peak);
    TF1* createGaussianFit(int maxBin);
    int findMaxBin();
    void detectAndFitPeaks();
    bool checkPredictedEnergies(double predictedEnergy, const double knownEnergies[], int size, float errorAdmitted) const;
    double refineCalibration(const double knownEnergies[], int size) const;
    void findStartOfPeak(Peak& peak);
    void initializeCalibratedHist();
    double getInterpolatedContent(int bin_original, double binCenter_original) const;
public:
    // Constructori și Destructor
    Histogram(int xMin, int xMax, int maxFWHM, int numberOfPeaks, TH1D *mainHist);
    ~Histogram();

    // Funcții publice
    void findPeaks();
    bool checkConditions(const Peak& peak) const;
    void outputPeaksDataJson(std::ofstream &file);
    void applyXCalibration();
    void calibratePeaks(const double knownEnergies[], int size);
    void printHistogramWithPeaksRoot(TFile *outputFile) const;
    void printCalibratedHistogramRoot(TFile *outputFile) const;
};

#endif // HISTOGRAM_H
