#include "../include/ArgumentsManager.h"
#include "../include/ErrorHandle.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <filesystem> 
#include <cctype>

ArgumentsManager::ArgumentsManager(int argc, char *argv[])
    : energyFilePath(getDataFolderPath()),
      energyProcessor(energyFilePath)
{
    parseArguments(argc, argv);
}

void ArgumentsManager::parseArguments(int argc, char *argv[])
{
    if (argc < 2)
    {
        printUsage();
        return;
    }

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-hf" || arg == "--histogram_file")
        {
            std::string input = argv[++i];
            if (isNumber(input))
            {
                int runNumber = std::stoi(input);
                std::string filename = getHistogramFilename(runNumber);
                if (!filename.empty())
                {
                    histogramFilePath = filename;
                }
                else
                {
                    ErrorHandle::getInstance().logStatus("Histogram file for run number " + std::to_string(runNumber) + " not found. Please provide a valid filename.");
                }
            }
            else
            {
                histogramFilePath = input;
            }
        }
        else if (arg == "-hn" || arg == "--histogram_name")
        {
            TH2histogram_name = argv[++i];
        }
        else if (arg == "-ef" || arg == "--energy_file")
        {
            // energyFilePath = argv[++i];
            // energyProcessor = sortEnergy(energyFilePath); // Update energyProcessor if the energy file changes
        }
        else if (arg == "-l" || arg == "-limits")
        {
            Xmin = std::stof(argv[++i]);
            Xmax = std::stof(argv[++i]);
            MinAmplitude = std::stof(argv[++i]);
            MaxAmplitude = std::stof(argv[++i]);
            FWHMmax = std::stof(argv[++i]);
        }
        else if (arg == "-sp" || arg == "--save_path")
        {
            savePath = argv[++i];
        }
        else if (arg == "-dt" || arg == "-detType")
        {
            detTypeStandard = std::stoi(argv[++i]);
        }
        else if (arg == "-se" || arg == "-serial")
        {
            serialStandard = argv[++i];
        }
        else if (arg == "-s" || arg == "-sources")
        {
            userInterfaceStatus = false;
            energyProcessor.chooseSources(i + 1, argc, argv);
            number_of_peaks = energyProcessor.getNumberOfPeaks();
            getSourcesNameRun();
            break;
        }
        else if (arg == "-j" || arg == "-json")
        {
            inputJsonFile = argv[++i];
            parseJsonFile();
        }
        else if (arg == "-d" || arg == "-domainLimits")
        {
            xMinDomain = std::stoi(argv[++i]);
            xMaxDomain = std::stoi(argv[++i]);
        }
        else if (arg == "-c" || arg == "-calib")
        {
            polynomialFitThreshold = std::stod(argv[++i]);
        }
        else if (arg == "-h" || arg == "--help")
        {
            printUsage();
            i = argc; // Stop parsing
            return;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << '\n';
            printUsage();
            return;
        }
    }
}

bool ArgumentsManager::parseNumericArgument(const char *arg, float &value, float min, float max)
{
    try
    {
        value = std::stof(arg);
        if (value >= min && value <= max)
        {
            return true;
        }
        std::cerr << "Value " << value << " out of range [" << min << ", " << max << "]\n";
    }
    catch (...)
    {
        std::cerr << "Failed to parse numeric value: " << arg << '\n';
    }
    return false;
}

// Add helper function to check if a string is a number
bool ArgumentsManager::isNumber(const std::string &s) const
{
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Modify function to search current directory for .root file containing run number
std::string ArgumentsManager::getHistogramFilename(int runNumber) const
{
    std::string runStr = std::to_string(runNumber);
    for (const auto &entry : std::filesystem::directory_iterator("."))
    {
        if (entry.is_regular_file())
        {
            std::string filename = entry.path().filename().string();
            if (filename.find(runStr) != std::string::npos && filename.find(".root") != std::string::npos)
            {
                return filename;
            }
        }
    }
    return "";
}

void ArgumentsManager::printUsage() const
{
    std::cout << "Usage: program [options]\n"
              << "Options:\n"
              << "  -hf, --histogram_file <path or run number>    Set histogram file path or run number\n" // Updated description
              << "  -hn, --histogram_name <name>                  Set histogram name\n"
              << "  -l, -limits <min> <max> <minA> <maxA> <fwhm>    Set analysis limits\n"
              << "  -sp, --save_path <path>                       Set save path\n"
              << "  -dt, -detType <type>                          Set detector type\n"
              << "  -s, -sources <source...>                      Specify sources\n"
              << "  -j, -json <file>                              Specify JSON configuration file\n"
              << "  -d, -domainLimits <min> <max>                  Set domain limits\n"
              << "  -c, -calib <threshold>                        Set calibration threshold\n";
}

std::string ArgumentsManager::getExecutableDir() const
{
    char buffer[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);

    if (len != -1)
    {
        buffer[len] = '\0'; // Null-terminate the path
        std::string fullPath = buffer;

        // Extract directory path by removing the executable name
        return fullPath.substr(0, fullPath.find_last_of('/'));
    }
    else
    {
        ErrorHandle::getInstance().logStatus("Failed to get executable directory.");
        return "";
    }
}

std::string ArgumentsManager::getDataFolderPath() const
{
    std::string exeDir = getExecutableDir();
    if (!exeDir.empty())
    {
        return exeDir + "/data/calibration_sources.json";
    }
    else
    {
        return "data/calibration_sources.json";
    }
}

bool ArgumentsManager::isDomainLimitsSet() const
{
    return xMinDomain != -1 && xMaxDomain != -1;
}

int ArgumentsManager::getNumberColumnSpecified(int histogramNumber) const
{
    auto it = std::find(domain.begin(), domain.end(), histogramNumber);
    if (it != domain.end())
    {
        return std::distance(domain.begin(), it);
    }
    else
    {
        return -1;
    }
}

bool ArgumentsManager::checkIfRunIsValid() const
{
    if (inputJsonFile.empty())
    {
        return false;
    }
    else
    {
        return true;
    }
}

void ArgumentsManager::getSourcesNameRun()
{
    std::string sourcesName;
    const auto &usedSources = energyProcessor.getRequestedSources();
    for (size_t i = 0; i < usedSources.size(); i++)
    {
        if (i != 0)
        {
            sourcesName += "+";
        }
        sourcesName += usedSources[i];
    }
    this->sourceName = sourcesName;
}

void ArgumentsManager::printAllArguments() const
{
    std::cout << "Number of peaks: " << number_of_peaks << std::endl;
    std::cout << "Histogram file path: " << histogramFilePath << std::endl;
    std::cout << "Histogram name: " << TH2histogram_name << std::endl;
    std::cout << "Energy file path: " << energyFilePath << std::endl;
    std::cout << "Xmin: " << Xmin << std::endl;
    std::cout << "Xmax: " << Xmax << std::endl;
    std::cout << "FWHMmax: " << FWHMmax << std::endl;
    std::cout << "Min amplitude: " << MinAmplitude << std::endl;
    std::cout << "Max amplitude: " << MaxAmplitude << std::endl;
    std::cout << "Save path: " << savePath << std::endl;
    std::cout << "Sources: " << getSourcesName() << std::endl;
}

void ArgumentsManager::setNumberOfPeaks(int peaks)
{
    number_of_peaks = peaks;
}

void ArgumentsManager::parseJsonFile()
{
    std::ifstream file(inputJsonFile);
    if (!file.is_open())
    {
        std::cerr << "Could not open Lut file: " << inputJsonFile << std::endl;
        return;
    }

    std::string line;
    bool inObject = false;
    // Temporary variables for current object
    int tempDomain = -1;
    int valuesRead = 0;
    int tempDetType = detTypeStandard;
    std::string tempSerial = serialStandard;
    int tempAmpl = MinAmplitude;
    int tempFwhm = FWHMmax;
    fitLimits tempLimits = {static_cast<int>(Xmin), static_cast<int>(Xmax)};
    PTLimits tempPTLimits = {static_cast<int>(MinAmplitude), static_cast<int>(MaxAmplitude)};

    while (std::getline(file, line))
    {
        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

        // Start of new object
        if (line.find("{") != std::string::npos)
        {
            inObject = true;
            // Reset temporary values to defaults
            tempDomain = -1;
            tempDetType = detTypeStandard;
            tempSerial = serialStandard;
            tempAmpl = MaxAmplitude;
            tempFwhm = FWHMmax;
            tempLimits = {static_cast<int>(Xmin), static_cast<int>(Xmax)};
            tempPTLimits = {static_cast<int>(MinAmplitude), static_cast<int>(MaxAmplitude)};
            continue;
        }

        // End of object - push all values to vectors
        if (line.find("}") != std::string::npos && inObject)
        {
            if (tempDomain != -1)
            { // Only add if domain was specified
                domain.push_back(tempDomain);
                detType.push_back(tempDetType);
                serial.push_back(tempSerial);
                ampl.push_back(tempAmpl);
                fwhm.push_back(tempFwhm);
                limits.push_back(tempLimits);
                ptLimits.push_back(tempPTLimits);
            }
            inObject = false;
            continue;
        }

        if (!inObject)
            continue;

        // Parse values within object
        if (line.find("\"domain\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                tempDomain = std::stoi(line.substr(colonPos + 1));
                valuesRead++;
            }
        }
        else if (line.find("\"detType\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                tempDetType = std::stoi(line.substr(colonPos + 1));
                valuesRead++;
            }
        }
        else if (line.find("\"serial\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                std::string value = line.substr(colonPos + 1);
                value.erase(std::remove(value.begin(), value.end(), '\"'), value.end());
                tempSerial = value;
                valuesRead++;
            }
        }
        else if (line.find("\"ampl\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                tempAmpl = std::stoi(line.substr(colonPos + 1));
                valuesRead++;
            }
        }
        else if (line.find("\"fwhm\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                tempFwhm = std::stoi(line.substr(colonPos + 1));
                valuesRead++;
            }
        }
        else if (line.find("\"fitLimits\"") != std::string::npos)
        {
            std::getline(file, line);
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            tempLimits.Xmin = std::stoi(line.substr(line.find("[") + 1));

            std::getline(file, line);
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            tempLimits.Xmax = std::stoi(line.substr(0, line.find("]")));
            valuesRead += 2;
        }
        else if (line.find("\"PTLimits\"") != std::string::npos)
        {
            std::getline(file, line);
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            tempPTLimits.MinAmplitude = std::stoi(line.substr(line.find("[") + 1));

            std::getline(file, line);
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            tempPTLimits.MaxAmplitude = std::stoi(line.substr(0, line.find("]")));
            valuesRead += 2;
        }
    }
    ErrorHandle::getInstance().logStatus("Data read successfully from Lut. Data read: " + std::to_string(valuesRead));
    file.close();
}

void ArgumentsManager::printArgumentsInput() const
{
    for (int i = 0; i < domain.size(); i++)
    {
        std::cout << "Domain: " << domain[i] << std::endl;
        std::cout << "DetType: " << detType[i] << std::endl;
        std::cout << "Serial: " << serial[i] << std::endl;
        std::cout << "Ampl: " << ampl[i] << std::endl;
        std::cout << "Fwhm: " << fwhm[i] << std::endl;
        std::cout << "FitLimits: " << limits[i].Xmin << " " << limits[i].Xmax << std::endl;
        std::cout << "PTLimits: " << ptLimits[i].MinAmplitude << " " << ptLimits[i].MaxAmplitude << std::endl;
    }
}

bool ArgumentsManager::fileExists(const std::string &path) const
{
    std::ifstream f(path);
    return f.good();
}
