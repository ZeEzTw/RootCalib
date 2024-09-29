#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>                 // pentru strcmp
#include "../include/sortEnergy.h" // Include sortEnergy header

class ArgumentsManager
{
private:
    int number_of_peaks = 1; // Default value
    std::string histogramFilePath = "data/data.root";
    std::string TH2histogram_name = "mDelila_raw";
    std::string energyFilePath = "data/calibration_sources.json";
    float Xmin = 0.0f;
    float Xmax = 1000000.0f;
    float FWHMmax = 1000.0f;
    float MinAmplitude = 0.0f;
    float MaxAmplitude = 1000000000.0f;
    std::string savePath = "output/";
    bool userInterfaceStatus = true; // Implicită pornirea UI-ului
    std::vector<std::string> userdSources;
    sortEnergy energyProcessor; // Changed to an object instead of a reference

public:
    ArgumentsManager(int argc, char *argv[]);

    void parseArguments(int argc, char *argv[]);
    std::string getSourcesName(); // Metodă care returnează numele surselor
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
    sortEnergy getEnergyProcessor() { return energyProcessor; }
    std::string getSavePath() const { return savePath; }
};
