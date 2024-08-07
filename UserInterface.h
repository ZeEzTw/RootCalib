#include <vector>
#include <string>
#include <iostream>
#include "Histogram.h"
#include "sortEnergy.h"

class UserInterface
{
public:
    void askAboutPeaks(std::vector<Histogram> &histograms, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated);
    double *askAboutSource(sortEnergy &energys, int &size);
    void showCalibrationInfo(const Histogram &histogram);
    };