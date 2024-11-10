#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>              
#include "../include/sortEnergy.h" 
#include <unistd.h>
#include <limits.h>

class ArgumentsManager
{
private:
    int number_of_peaks = 1;
    std::string histogramFilePath = "data/data.root";
    std::string TH2histogram_name = "mDelila_raw";
    std::string energyFilePath = "data/calibration_sources.json";
    float Xmin = 0.0f;
    float Xmax = 1000000.0f;
    float MinAmplitude = 0.0f;
    float MaxAmplitude = 1000000000.0f;
    float FWHMmax = 1000.0f;
    std::string savePath;
    int detTypeStandard = 2;
    std::string serialStandard = "CL";
    bool userInterfaceStatus = true; 
    std::vector<std::string> userdSources;
    sortEnergy energyProcessor;
    std::string sourceName;
    std::string inputJsonFile;
    int xMinDomain = -1;
    int xMaxDomain = -1;
    std::vector<int> domain;
    std::vector<int> detType;
    std::vector<std::string> serial;
    std::vector<int> ampl;
    std::vector<int> fwhm;
    float polynomialFitThreshold = 1e-3;

    struct fitLimits
    {
        int Xmin;
        int Xmax;
    };
    std::vector<fitLimits> limits;

    struct PTLimits
    {
        int MinAmplitude;
        int MaxAmplitude;
    };
    std::vector<PTLimits> ptLimits;

public:
    ArgumentsManager(int argc, char *argv[]);
    bool isDomainLimitsSet();
    int getNumberColumnSpecified(int histogramNumber);
    bool checkIfRunIsValid();
    void parseArguments(int argc, char *argv[]);
    void getSourcesNameRun(); 
    std::string getSourcesName() const { return sourceName; }
    void setSourceName(std::string &sourceName) { this->sourceName = sourceName; }
    bool isUserInterfaceEnabled() const { return userInterfaceStatus; }
    void printAllArguments();
    // Getters
    int getNumberOfPeaks() const { return number_of_peaks; }
    std::string getHistogramFilePath() const { return histogramFilePath; }
    std::string getHistogramName() const { return TH2histogram_name; }
    std::string getEnergyFilePath() const { return energyFilePath; }
    float getXmin() const { return Xmin; }
    float getXmax() const { return Xmax; }
    float getFWHMmax() const { return FWHMmax; }
    float getMinAmplitude() const { return MinAmplitude; }
    float getMaxAmplitude() const { return MaxAmplitude; }
    int getDetTypeStandard() const { return detTypeStandard; };
    std::string getSerialStandard() const { return serialStandard; };
    void printArgumentsInput();
    sortEnergy getEnergyProcessor() { return energyProcessor; }
    std::string getSavePath() const { return savePath; }
    void setNumberOfPeaks(int peaks);
    void parseJsonFile();
    int getXmaxDomain();
    int getXminDomain();
    int getXminFile(int position);
    int getXmaxFile(int position);
    int getFWHMmaxFile(int position);
    float getMaxAmplitudeFile(int position) const;
    std::string getSerialFile(int position) const;
    int getDetTypeFile(int position) const;
    std::string getHistogramNameFile(int position);
    float getPolynomialFitThreshold() const { return polynomialFitThreshold; }
    std::string getExecutableDir();
    std::string getDataFolderPath();
};
