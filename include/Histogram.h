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
#include <TLatex.h>
#include <string>
#include <TGraph.h>
class Histogram
{
private:
    int xMin, xMax;
    int maxFWHM;
    float minAmplitude;
    float maxAmplitude;
    int numberOfPeaks;
    TH1D *mainHist;
    TH1D *tempHist;
    TH1D *calibratedHist;
    float m, b;
    int polinomDegree;
    unsigned int peakMatchCount;
    std::vector<Peak> peaks;
    int peakCount;
    std::string TH2histogram_name;
    std::string sourceName;
    // Funcții private
    void eliminatePeak(const Peak &peak);
    TF1 *createGaussianFit(int maxBin);
    int findMaxBin();
    void detectAndFitPeaks();
    bool checkPredictedEnergies(double predictedEnergy, const double knownEnergies[], int size, float errorAdmitted) const;
    double refineCalibration(const double knownEnergies[], int size) const;
    void findStartOfPeak(Peak &peak);
    void initializeCalibratedHist();
    double getInterpolatedContent(int bin_original, double binCenter_original) const;
    int getTheDegreeOfPolynomial() const;

public:
    // Constructori și Destructor
    Histogram(int xMin, int xMax, int maxFWHM, float minAmplitude, float maxAmplitude, int numberOfPeaks, TH1D *mainHist, const std::string &TH2histogram_name, std::string sourceName);
    Histogram(const Histogram &histogram);
    ~Histogram();
    Histogram &operator=(const Histogram &histogram);
    // Funcții publice
    void findPeaks();
    bool checkConditions(const Peak &peak) const;
    void outputPeaksDataJson(std::ofstream &file);
    void applyXCalibration();
    void calibratePeaks(const double knownEnergies[], int size);
    void changePeak(int peakNumber, double newPosition);
    void printHistogramWithPeaksRoot(TFile *outputFile);
    void printCalibratedHistogramRoot(TFile *outputFile) const;
    const char *returnNameOfHistogram() const;
    unsigned int getpeakMatchCount() const;
    TH1D *getCalibratedHist() const;
    TH1D *getMainHist() const;
};

#endif // HISTOGRAM_H
