#ifndef ERRORHANDLE_H
#define ERRORHANDLE_H

#include <string>
#include <vector>

// Structure to represent an error entry
struct ErrorEntry
{
    std::string timestamp;
    int error;
    std::string error_message;
    std::string error_solution;
};

// Structure to represent a status update
struct StatusEntry
{
    std::string timestamp;
    std::string message;
};

class ErrorHandle
{
public:
    enum ErrorCode
    {
        SUCCESS = 0,
        TOO_FEW_ARGUMENTS = 1,
        INVALID_SOURCE_NAMES = 2,
        INVALID_INPUT_FILE = 3,
        INVALID_TH2F_HISTOGRAM = 4,
        NO_VALID_PEAKS = 5,           // Added missing comma
        INVALID_OUTPUT_FILE = 6,      // Added missing comma
        NO_PEAKS_FOR_CALIBRATION = 7, // Added missing comma
        LUT_FILE_NOT_FOUND = 8        // Last enumerator can optionally have a comma
    };

    // Static method to get the Singleton instance
    static ErrorHandle &getInstance();

    // Delete copy constructor and assignment operator
    ErrorHandle(const ErrorHandle &) = delete;
    ErrorHandle &operator=(const ErrorHandle &) = delete;

    // Member function declarations
    void errorHandle(int errorNumber);
    void saveLogFile();
    void logStatus(const std::string &statusMessage);
    void logLutFileInput(const std::string &lutFileName, int rowsRead);
    void logArrayWithCalibratedValues(const double *array, int size);
    void startProgram();

    // Setter methods for configuration
    void setUserInterfaceActive(bool isActive);
    void setPathForSave(const std::string &path);

private:
    // Private constructor
    ErrorHandle();

    bool isUserInterfaceActive;
    std::string pathForSave;
    std::vector<ErrorEntry> errors;
    std::vector<StatusEntry> status_updates;
    bool isInputFileValid;
    bool isDataReadSuccessful;

    // Helper function declarations
    std::string getCurrentTime();
    void writeProblemToJsonErrorFile(int errorNumber, const std::string &errorMessage, const std::string &errorSolution);
};

#endif // ERRORHANDLE_H