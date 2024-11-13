/**
 * @class ArgumentsManager
 * @brief Handles parsing and management of command-line arguments and calibration configurations
 *
 * The ArgumentsManager class is responsible for:
 * - Parsing and validating command-line arguments for the calibration process
 * - Managing calibration data parameters, including limits and thresholds
 * - Processing JSON configuration files to set up and manage calibration runs
 * - Handling detector configurations and managing energy calibration data
 * - Controlling the activation status of the user interface
 * - Generating calibrated energy arrays for the calibration process
 * - Providing validated and organized data to the TaskHandler for further processing
 */
#include <string>
#include "../include/CalibrationDataProvider.h"

class ArgumentsManager
{
private:
    // File data structures
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

    // Core configuration
    std::string histogramFilePath = "data/data.root";
    std::string TH2histogram_name = "mDelila_raw";
    std::string energyFilePath = "data/calibration_sources.json";
    std::string savePath;
    std::string inputJsonFile;
    std::string sourceName;
    CalibrationDataProvider energyProcessor;

    // Analysis parameters
    int Xmin = 0;
    int Xmax = 1000000;
    float MinAmplitude = 0.0f;
    float MaxAmplitude = 1e6f;
    float FWHMmax = 1e4f;
    float polynomialFitThreshold = 1e-6f;
    int number_of_peaks = 1;

    // Domain configuration
    int xMinDomain = -1;
    int xMaxDomain = -1;
    std::vector<int> domain;

    // Detector configuration
    int detTypeStandard = 2;
    std::string serialStandard = "CL";
    std::vector<int> detType;
    std::vector<std::string> serial;

    // Analysis vectors
    std::vector<int> ampl;
    std::vector<int> fwhm;
    std::vector<fitLimits> limits;
    std::vector<PTLimits> ptLimits;
    std::vector<std::string> usedSources;

    // State
    bool userInterfaceStatus = true;

    // Private helper methods
    bool validateInputParameters() const;
    bool parseNumericArgument(const char *arg, float &value, float min, float max);
    std::string getHistogramFilename(int runNumber) const;
    bool isNumber(const std::string &s) const;
    bool fileExists(const std::string &path) const;

public:
    // Constructor and main interface
    explicit ArgumentsManager(int argc, char *argv[]);

    // Main operations
    void parseArguments(int argc, char *argv[]);
    void parseJsonFile();
    void getSourcesNameRun();

    // Validation and status
    bool isDomainLimitsSet() const;
    bool checkIfRunIsValid() const;
    bool isUserInterfaceEnabled() const { return userInterfaceStatus; }

    // Print functions
    void printUsage() const;
    void printAllArguments() const;
    void printArgumentsInput() const;

    // Getters for source information
    std::string getSourcesName() const { return sourceName; }
    void setSourceName(const std::string &name) { sourceName = name; }

    // Getters and utilities
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

    // Domain and file getters
    int getXmaxDomain() const { return xMaxDomain; }
    int getXminDomain() const { return xMinDomain; }
    int getXminFile(int position) const { return limits[position].Xmin; }
    int getXmaxFile(int position) const { return limits[position].Xmax; }
    int getFWHMmaxFile(int position) const { return fwhm[position]; }
    int getMinAmplitudeFile(int position) const { return ampl[position]; }
    std::string getSerialFile(int position) const { return serial[position]; }
    int getDetTypeFile(int position) const { return detType[position]; }
    std::string getHistogramNameFile(int position) const { return serial[position]; }

    // Declare non-inlined functions
    int getNumberColumnSpecified(int histogramNumber) const;

    // Other functions
    void setNumberOfPeaks(int peaks);
    std::string getExecutableDir() const;
    std::string getDataFolderPath() const;
};
