/**
 * @brief Manages command-line arguments and calibration configurations for the Root Calibration system
 *
 * This file contains the ArgumentsManager class which handles all configuration aspects
 * of the calibration process, including command-line parsing (main), file management(lut file),
 * and parameter validation(used all over the program).
 */

#include <string>
#include "../include/CalibrationDataProvider.h"

class ArgumentsManager
{
private:
    //--------------------
    // Limits and Ranges
    //--------------------
    struct fitLimits
    {
        int Xmin;
        int Xmax;
    };

    struct PTLimits
    {
        int MinAmplitude;
        int MaxAmplitude;
    };

    //--------------------
    // File Configuration and locations
    //--------------------
    std::string histogramFilePath = "data/data.root";
    std::string TH2histogram_name = "mDelila_raw";
    std::string energyFilePath = "data/calibration_sources.json";
    std::string savePath;
    std::string inputJsonFile;
    std::string sourceName;
    CalibrationDataProvider energyProcessor;

    //--------------------
    // Analysis Parameters
    //--------------------
    int Xmin = 0;
    int Xmax = 1000000;
    float MinAmplitude = 0.0f;
    float MaxAmplitude = 1e10f;
    float FWHMmax = 1e4f;
    float polynomialFitThreshold = 1e-6f;
    int number_of_peaks = 1;

    //--------------------
    // Domain Configuration
    //--------------------
    int xMinDomain = -1;
    int xMaxDomain = -1;
    std::vector<int> domain;

    //--------------------
    // Detector Configuration
    //--------------------
    int detTypeStandard = 2;
    std::string serialStandard = "CL";
    std::vector<int> detType;
    std::vector<std::string> serial;

    //--------------------
    // Analysis Vectors
    //--------------------
    std::vector<int> ampl;
    std::vector<int> fwhm;
    std::vector<fitLimits> limits;
    std::vector<PTLimits> ptLimits;
    std::vector<std::string> usedSources;

    //--------------------
    // State if user interface is on
    //--------------------
    bool userInterfaceStatus = true;

    //--------------------
    // Private Helper Methods
    //--------------------
    bool validateInputParameters() const;
    bool parseNumericArgument(const char *arg, float &value, float min, float max);
    std::string getHistogramFilename(int runNumber) const;
    bool isNumber(const std::string &s) const;
    bool fileExists(const std::string &path) const;

public:
    //--------------------
    // Constructor and Main Interface
    //--------------------
    explicit ArgumentsManager(int argc, char *argv[]);

    //--------------------
    // Main Operations
    //--------------------
    void parseArguments(int argc, char *argv[]);
    void parseJsonFile();
    void getSourcesNameRun();

    //--------------------
    // Validation and Status
    //--------------------
    bool isDomainLimitsSet() const;
    bool checkIfRunIsValid() const;
    bool isUserInterfaceEnabled() const { return userInterfaceStatus; }

    //--------------------
    // Print Functions
    //--------------------
    void printUsage() const;
    void printAllArguments() const;
    void printArgumentsInput() const;

    //--------------------
    // Getters for Source Information
    //--------------------
    std::string getSourcesName() const { return sourceName; }
    void setSourceName(const std::string &name) { sourceName = name; }

    //--------------------
    // Getters and Utilities
    //--------------------
    int getNumberOfPeaks() const { return number_of_peaks; }
    std::string getHistogramFilePath() const { return histogramFilePath; }
    std::string getHistogramName() const { return TH2histogram_name; }
    std::string getEnergyFilePath() const { return energyFilePath; }
    int getXmin() const { return Xmin; }
    int getXmax() const { return Xmax; }
    float getFWHMmax() const { return FWHMmax; }
    float getMinAmplitude() const { return MinAmplitude; }
    float getMaxAmplitude() const { return MaxAmplitude; }
    int getDetTypeStandard() const { return detTypeStandard; }
    std::string getSerialStandard() const { return serialStandard; }
    CalibrationDataProvider getEnergyProcessor() const { return energyProcessor; }
    std::string getSavePath() const { return savePath; }
    float getPolynomialFitThreshold() const { return polynomialFitThreshold; }

    //--------------------
    // Domain and File Getters
    //--------------------
    int getXmaxDomain() const { return xMaxDomain; }
    int getXminDomain() const { return xMinDomain; }
    int getXminFile(int position) const { return limits[position].Xmin; }
    int getXmaxFile(int position) const { return limits[position].Xmax; }
    int getFWHMmaxFile(int position) const { return fwhm[position]; }
    int getMinAmplitudeFile(int position) const { return ampl[position]; }
    std::string getSerialFile(int position) const { return serial[position]; }
    int getDetTypeFile(int position) const { return detType[position]; }
    std::string getHistogramNameFile(int position) const { return serial[position]; }

    //--------------------
    // Declare Non-Inlined Functions
    //--------------------
    int getNumberColumnSpecified(int histogramNumber) const;

    //--------------------
    // Other Functions
    //--------------------
    void setNumberOfPeaks(int peaks);
    std::string getExecutableDir() const;
    std::string getDataFolderPath() const;
};
