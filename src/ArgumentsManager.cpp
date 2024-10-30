#include "../include/ArgumentsManager.h"
#include <iostream>
#include <cstring> // pentru strcmp

ArgumentsManager::ArgumentsManager(int argc, char *argv[])
    : energyProcessor(energyFilePath)
{ // Initialize energyProcessor with the energyFilePath
    parseArguments(argc, argv);
}

void ArgumentsManager::parseArguments(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-hf" || arg == "--histogram_file")
        {
            histogramFilePath = argv[++i];
        }
        else if (arg == "-hn" || arg == "--histogram_name")
        {
            TH2histogram_name = argv[++i];
        }
        else if (arg == "-ef" || arg == "--energy_file")
        {
            energyFilePath = argv[++i];
            energyProcessor = sortEnergy(energyFilePath); // Update energyProcessor if the energy file changes
        }
        else if (arg == "-limits")
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
        else if (arg == "-detType")
        {
            detTypeStandard = std::stoi(argv[++i]);
        }
        else if (arg == "-serial")
        {
            serialStandard = argv[++i];
        }
        else if (arg == "-sources")
        {
            // std::cout<<"Aici"<<std::endl;
            userInterfaceStatus = false; // Dacă se specifică sursele, UI-ul se oprește
            energyProcessor.chooseSources(i + 1, argc, argv);
            number_of_peaks = energyProcessor.getNumberOfPeaks();
            // std::cout<<"start: "<<i<<std::endl;
            // std::cout<<"argc: "<<argc<<std::endl;
            // std::cout<<"Sources: "<<argv[i + 1]<<std::endl;
            break; // Nu mai continuăm, deoarece sursele sunt ultimele argumente
        }
        else if (arg == "-json")
        {
            inputJsonFile = argv[++i];
            // parseJsonFile();
        }
        else if (arg == "-domainLimits")
        {
            xMinDomain = std::stoi(argv[++i]);
            xMaxDomain = std::stoi(argv[++i]);
        }
        else if (arg == "-calib")
        {
            polynomialFitThreshold = std::stod(argv[++i]);
        }
    }
}

bool ArgumentsManager::isDomainLimitsSet()
{
    return xMinDomain != -1 && xMaxDomain != -1;
}

int ArgumentsManager::GetNumberColomSpecified(int histogramNumber)
{
    auto it = std::find(domain.begin(), domain.end(), histogramNumber);

    if (it != domain.end())
    {
        // std::cout << "Found at position: " << std::distance(domain.begin(), it) << std::endl;
        return std::distance(domain.begin(), it);
    }
    else
    {
        return -1;
    }
}

bool ArgumentsManager::checkIfRunIsValid()
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
    const auto &usedSources = energyProcessor.getRequestedSources(); // Accesează sursele corect
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

void ArgumentsManager::printAllArguments()
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
        // std::cerr << "Could not open file: " << inputJsonFile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line))
    {
        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

        if (line.find("\"domain\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                domain.push_back(std::stoi(line.substr(colonPos + 1)));
            }
        }

        if (line.find("\"detType\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                detType.push_back(std::stoi(line.substr(colonPos + 1)));
            }
        }

        if (line.find("\"serial\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                std::string value = line.substr(colonPos + 1);
                value.erase(remove(value.begin(), value.end(), '\"'), value.end()); // Îndepărtăm ghilimelele
                serial.push_back(value);
            }
        }

        if (line.find("\"ampl\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                ampl.push_back(std::stoi(line.substr(colonPos + 1)));
            }
        }

        if (line.find("\"fwhm\"") != std::string::npos)
        {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos)
            {
                fwhm.push_back(std::stoi(line.substr(colonPos + 1)));
            }
        }

        if (line.find("\"fitLimits\"") != std::string::npos)
        {
            fitLimits temp;
            std::getline(file, line); // Citim linia pentru primul element (Xmin)
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            temp.Xmin = std::stoi(line.substr(line.find("[") + 1)); // Valoarea Xmin

            std::getline(file, line); // Citim linia pentru al doilea element (Xmax)
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            temp.Xmax = std::stoi(line.substr(0, line.find("]"))); // Valoarea Xmax

            limits.push_back(temp);
        }

        // Parse 'PTLimits' (listă de două valori int)
        if (line.find("\"PTLimits\"") != std::string::npos)
        {
            PTLimits temp;
            std::getline(file, line); // Citim linia pentru MinAmplitude
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            temp.MinAmplitude = std::stoi(line.substr(line.find("[") + 1)); // Valoarea MinAmplitude

            std::getline(file, line); // Citim linia pentru MaxAmplitude
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
            temp.MaxAmplitude = std::stoi(line.substr(0, line.find("]"))); // Valoarea MaxAmplitude

            ptLimits.push_back(temp);
        }
    }

    file.close();
}

void ArgumentsManager::printArgumentsInput()
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

int ArgumentsManager::getXmaxDomain()
{
    return xMaxDomain;
}

int ArgumentsManager::getXminDomain()
{
    return xMinDomain;
}

int ArgumentsManager::getXminFile(int position)
{
    return limits[position].Xmin;
}

int ArgumentsManager::getXmaxFile(int position)
{
    return limits[position].Xmax;
}

int ArgumentsManager::getFWHMmaxFile(int position)
{
    return fwhm[position];
}

float ArgumentsManager::getMaxAmplitudeFile(int position) const
{
    return ampl[position];
}

std::string ArgumentsManager::getHistogramNameFile(int position)
{
    return serial[position];
}

int ArgumentsManager::getDetTypeFile(int position) const
{
    return detType[position];
}

std::string ArgumentsManager::getSerialFile(int position) const
{
    return serial[position];
}
