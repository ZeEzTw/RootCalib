#include "../include/ErrorHandle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring> // For strerror

// Singleton instance retrieval
ErrorHandle &ErrorHandle::getInstance()
{
    static ErrorHandle instance;
    return instance;
}

// Private constructor implementation
ErrorHandle::ErrorHandle()
    : isUserInterfaceActive(false), pathForSave(""), isInputFileValid(false), isDataReadSuccessful(false)
{
    // Initialize containers if necessary
}

// Setter for isUserInterfaceActive
void ErrorHandle::setUserInterfaceActive(bool isActive)
{
    isUserInterfaceActive = isActive;
}

// Setter for pathForSave
void ErrorHandle::setPathForSave(const std::string &path)
{
    pathForSave = path;
}

// Function to get current time as string in HH:MM:SS format
std::string ErrorHandle::getCurrentTime()
{
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%H:%M:%S");
    return ss.str();
}

void ErrorHandle::errorHandle(int errorNumber)
{
    std::string errorMessage;
    std::string errorSolution;

    switch (errorNumber)
    {
    case SUCCESS:
        errorMessage = "Program finished successfully.";
        errorSolution = "";
        break;
    case TOO_FEW_ARGUMENTS:
        errorMessage = "Too few arguments to run.";
        errorSolution = "Please provide the necessary arguments. user -h for help.";
        break;
    case INVALID_SOURCE_NAMES:
        errorMessage = "Source names for calibration are not valid.";
        errorSolution = "Check the spelling of the source names or open then file calibration_sources.json in data to see spelling or run the program witout soruces name and UserInterface will show you the available sources.";
        break;
    case INVALID_INPUT_FILE:
        errorMessage = "Input file is invalid or cannot be opened.";
        errorSolution = "Check the spelling of the input file path.";
        break;
    case INVALID_OUTPUT_FILE:
        errorMessage = "Output file is invalid or cannot be opened.";
        errorSolution = "Check the spelling of the output file path, if this file is not specify, it will be created automatically. "
                        "Fired by: Could not create output file for histograms or Could not create output file for calibrated histograms or Could not create output file for combined histograms.";
        break;
    case INVALID_TH2F_HISTOGRAM:
        errorMessage = "The TH2F histogram with data is not valid.";
        errorSolution = "Check the spelling of the histogram name or Input File path is invalid, if not specified it will be m_delida.";
        break;
    case NO_VALID_PEAKS:
        errorMessage = "No valid peaks for calibration (number of peaks = 0).";
        errorSolution = "Ensure that there are valid peaks in the data.";
        break;
    case NO_PEAKS_FOR_CALIBRATION:
        errorMessage = "No peaks available for calibration. n = 0";
        errorSolution = "Probably no peaks pass the conditions (check LUT file, input data from terminal or default conditions).";
        break;
    case LUT_FILE_NOT_FOUND:
        errorMessage = "LUT file not found. Called in ArgumentsManager::parseJsonFile()";
        errorSolution = "Check the spelling of the LUT file path.";
        break;
    default:
        errorMessage = "Unknown error occurred.";
        errorSolution = "Please consult the documentation.";
        break;
    }

    if (isUserInterfaceActive)
    {
        std::cerr << "Error Code: " << errorNumber << std::endl;
        std::cerr << "ErrorMessage: " << errorMessage << std::endl;
        if (!errorSolution.empty())
        {
            std::cerr << "Solution: " << errorSolution << std::endl;
        }
    }
    else
    {
        std::cerr << "Error Code " << errorNumber << std::endl;
    }

    // Log error to internal container with timestamp
    writeProblemToJsonErrorFile(errorNumber, errorMessage, errorSolution);
}

void ErrorHandle::writeProblemToJsonErrorFile(int errorNumber, const std::string &errorMessage, const std::string &errorSolution)
{
    ErrorEntry entry;
    entry.timestamp = getCurrentTime();
    entry.error = errorNumber;
    entry.error_message = errorMessage;
    entry.error_solution = errorSolution;
    errors.push_back(entry);
}

void ErrorHandle::saveLogFile()
{
    std::cout<<"Saving log file"<<std::endl;
    // Ensure the save directory exists; create it if it doesn't
    if (pathForSave.empty())
    {
        pathForSave = ".";
        // Optionally log the creation of the directory
        logStatus("Using current directory for saving logs, pathForSave was invalid: " + pathForSave);
    }
    else
    {
        struct stat info;
        if (stat(pathForSave.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR))
        {
            if (mkdir(pathForSave.c_str(), 0777) != 0)
            {
                std::cerr << "Error creating directory: " << strerror(errno) << std::endl;
                logStatus(std::string("Failed to create log directory: ") + strerror(errno));
                return; // Exit the function if directory creation fails
            }
            else
            {
                // Optionally log the creation of the directory
                logStatus("Created log directory: " + pathForSave);
            }
        }
        std::cout<<"pathForSave: "<<pathForSave<<std::endl;
    }

    std::ofstream logFile(pathForSave + "/error_log.json");
    if (logFile.is_open())
    {
        logFile << "{\n";

        // Write errors
        logFile << "  \"errors\": [\n";
        for (size_t i = 0; i < errors.size(); ++i)
        {
            logFile << "    {\n";
            logFile << "      \"error code\": " << errors[i].error << ",\n";
            logFile << "      \"error_message\": \"" << errors[i].error_message << "\",\n";
            logFile << "      \"error_solution\": \"" << errors[i].error_solution << "\",\n";
            logFile << "      \"timestamp\": \"" << errors[i].timestamp << "\"\n";

            logFile << "    }";
            if (i != errors.size() - 1)
                logFile << ",";
            logFile << "\n";
        }
        logFile << "  ],\n";

        // Write status updates
        logFile << "  \"status_updates\": [\n";
        for (size_t i = 0; i < status_updates.size(); ++i)
        {
            logFile << "    {\n";
            logFile << "      \"message\": \"" << status_updates[i].message << "\",\n";
            logFile << "      \"timestamp\": \"" << status_updates[i].timestamp << "\"\n";

            logFile << "    }";
            if (i != status_updates.size() - 1)
                logFile << ",";
            logFile << "\n";
        }
        logFile << "  ]\n";

        logFile << "}\n";
        logFile.close();

        // Optionally log that the file was saved successfully
        logStatus("Log file saved successfully at: " + pathForSave + "/error_log.json");
    }
    else
    {
        std::cerr << "Failed to open log file for writing." << std::endl;
        logStatus("Failed to open log file for writing at: " + pathForSave + "/error_log.json");
    }
}

void ErrorHandle::logStatus(const std::string &statusMessage)
{
    if (isUserInterfaceActive)
    {
        std::cout << statusMessage << std::endl;
    }
    // Log status to internal container with timestamp
    StatusEntry entry;
    entry.timestamp = getCurrentTime();
    entry.message = statusMessage;
    status_updates.push_back(entry);
}

void ErrorHandle::startProgram()
{
    if (isUserInterfaceActive)
        std::cout << "Program started successfully." << std::endl;
    logStatus("Program started successfully.");
}
void ErrorHandle::logLutFileInput(const std::string &lutFileName, int rowsRead)
{
    std::stringstream ss;
    ss << "The LUT file \"" << lutFileName << "\" was opened successfully.";
    ss << " " << rowsRead << " rows read from it.";
    logStatus(ss.str());
}

void ErrorHandle::logArrayWithCalibratedValues(const double *array, int size)
{
    std::stringstream ss;
    ss << "Calibrated energy array: [";
    for (int i = 0; i < size; ++i)
    {
        ss << array[i];
        if (i != size - 1)
            ss << ", ";
    }
    ss << "]";
    logStatus(ss.str());
}
