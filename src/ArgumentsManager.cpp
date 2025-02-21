#include "../include/ArgumentsManager.h"
#include "../include/ErrorHandle.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <dirent.h> // Include dirent.h for directory iteration
#include <cctype>
#include <nlohmann/json.hpp> // Include nlohmann/json

ArgumentsManager::ArgumentsManager(int argc, char *argv[])
    : energyFilePath(getDataFolderPath()),
      energyProcessor(energyFilePath)
{
    parseArguments(argc, argv); // take the data from the terminal when the class is initialized
}

void ArgumentsManager::parseArguments(int argc, char *argv[])
{
    if (argc < 2)//100% required -f / --histogram_file and -rp / --json: Path to the LUT JSON
    {
        printUsage();
        return;
    }
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-f" || arg == "--histogram_file")//location of the TH2F histogram file
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
        else if (arg == "-hn" || arg == "--histogram_name")//ex m_deliaRaw, name of the TH2F histogram
        {
            TH2histogram_name = argv[++i];
        }
        else if (arg == "-ef" || arg == "--energy_file")
        {
            // energyFilePath = argv[++i];
            // energyProcessor = sortEnergy(energyFilePath); // Update energyProcessor if the energy file changes
        }
        else if ((arg == "-l" || arg == "-limits") && i + 5 < argc)
        {
            Xmin = std::stoi(argv[++i]);
            Xmax = std::stoi(argv[++i]);
            MinAmplitude = std::stof(argv[++i]);
            MaxAmplitude = std::stof(argv[++i]);//not used really
            FWHMmax = std::stof(argv[++i]);
        }
        else if (arg == "-sp" || arg == "--save_path")//Output directory. Default: output/.
        {
            savePath = argv[++i];
        }
        else if (arg == "-dt" || arg == "-detType")//Detector type. Default: 2. Its updated by pyCalib
        {
            detTypeStandard = std::stoi(argv[++i]);
        }
        else if (arg == "-se" || arg == "-serial")//Detector serial number. Default: CL.
        {
            serialStandard = argv[++i];
        }
        //if not specified (here can pe puted mode then 1 soruces), the user interface will pop up
        //and the soruce need to be the last in the comand line
        else if (arg == "-s" || arg == "-sources")
        {
            userInterfaceStatus = false;
            energyProcessor.chooseSources(i + 1, argc, argv);
            number_of_peaks = energyProcessor.getNumberOfPeaks();
            getSourcesNameRun();
            break;
        }
        else if (arg == "-rp" || arg == "-json")//Path to the LUT JSON file
        {
            inputJsonFile = argv[++i];
        }
        else if (arg == "-sc" || arg == "-domainLimitsStart")//Domain limits for the TH2F histogram, 
        {
            xMinDomain = std::stoi(argv[++i]);
        }
        else if (arg == "-ec" || arg == "-domainLimitsEnd")
        {
            xMaxDomain = std::stoi(argv[++i]);
        }
        else if (arg == "-c" || arg == "-calib")//calibration polynomial threshold. Default: 1e-10.
        {
            polynomialFitThreshold = std::stod(argv[++i]);
        }
        else if (arg == "-h" || arg == "--help")//the help function, it will pop up when called or when not all requred data are given(ex: source)
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
    if (inputJsonFile.empty())
    {
        ErrorHandle::getInstance().logStatus("No JSON file provided. Please provide a JSON file.");
    }
    else
    {
        parseJsonFile();
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

// function to check if a string is a number, for only run number specified
bool ArgumentsManager::isNumber(const std::string &s) const
{
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// function to search current directory for .root file containing run number
std::string ArgumentsManager::getHistogramFilename(int runNumber) const
{
    std::string runStr = "_" + std::to_string(runNumber) + "_";
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(".")) != NULL)
    {
        while ((ent = readdir(dir)) != NULL)
        {
            std::string filename = ent->d_name;
            if (filename.find(runStr) != std::string::npos && filename.find(".root") != std::string::npos)
            {
                closedir(dir);
                return filename;
            }
        }
        closedir(dir);
    }
    return "";
}

//help function
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

//retrieves the directory path where the currently running executable is located
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
//This function determines the path where the calibration will be saved
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

//check if the run number is in the domain (in the lut)
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

//for the json, in case multiple sources are used, will handle them
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


/**
 * @brief Parses configuration data from a JSON file.
 * This method reads and processes a JSON file containing configuration data for detector parameters.
 * The method populates class member vectors with the parsed values.
 * Only entries with valid domain values (!= -1) are stored.
 * 
 * @throws May throw exceptions from file operations or JSON parsing
 * @note Logs success message upon completion via ErrorHandle singleton
 * rediding the data from the LUT JSON file
 */
void ArgumentsManager::parseJsonFile()
{
    std::ifstream file(inputJsonFile);
    if (!file.is_open())
    {
        std::cerr << "Could not open JSON file: " << inputJsonFile << std::endl;
        return;
    }

    nlohmann::json jsonData;
    file >> jsonData;

    for (const auto &item : jsonData)
    {
        int tempDomain = item.contains("domain") ? item["domain"].get<int>() : -1;
        int tempDetType = item.contains("detType") ? item["detType"].get<int>() : detTypeStandard;
        std::string tempSerial = item.contains("serial") ? item["serial"].get<std::string>() : serialStandard;
        int tempAmpl = item.contains("ampl") ? item["ampl"].get<int>() : MinAmplitude;
        int tempFwhm = item.contains("fwhm") ? item["fwhm"].get<int>() : FWHMmax;
        fitLimits tempLimits = {Xmin, Xmax};
        PTLimits tempPTLimits = {static_cast<int>(MinAmplitude), static_cast<int>(MaxAmplitude)};

        if (item.contains("fitLimits"))
        {
            tempLimits.Xmin = item["fitLimits"].contains("Xmin") ? item["fitLimits"]["Xmin"].get<int>() : Xmin;
            tempLimits.Xmax = item["fitLimits"].contains("Xmax") ? item["fitLimits"]["Xmax"].get<int>() : Xmax;
        }

        if (item.contains("PTLimits"))
        {
            tempPTLimits.MinAmplitude = item["PTLimits"].contains("MinAmplitude") ? item["PTLimits"]["MinAmplitude"].get<int>() : static_cast<int>(MinAmplitude);
            tempPTLimits.MaxAmplitude = item["PTLimits"].contains("MaxAmplitude") ? item["PTLimits"]["MaxAmplitude"].get<int>() : static_cast<int>(MaxAmplitude);
        }

        if (tempDomain != -1)
        {
            domain.push_back(tempDomain);
            detType.push_back(tempDetType);
            serial.push_back(tempSerial);
            ampl.push_back(tempAmpl);
            fwhm.push_back(tempFwhm);
            limits.push_back(tempLimits);
            ptLimits.push_back(tempPTLimits);
        }
    }

    ErrorHandle::getInstance().logStatus("Data read successfully from JSON file.");
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
