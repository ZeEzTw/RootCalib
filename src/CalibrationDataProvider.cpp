#include "../include/CalibrationDataProvider.h"
#include "../include/ErrorHandle.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
/**
 * @param filename Path to calibration configuration file
 */
CalibrationDataProvider::CalibrationDataProvider(const std::string &filename)
{
    ErrorHandle::getInstance().logStatus("Reading calibration data from: " + filename);
    parseJsonFile(filename);
}

CalibrationDataProvider &CalibrationDataProvider::operator=(const CalibrationDataProvider &other)
{
    if (this != &other)
    {
        sources = other.sources;
        energyMatrix = other.energyMatrix;
        requestedSources = other.requestedSources;
        numberOfPeaks = other.numberOfPeaks;
        probabilityMatrix = other.probabilityMatrix;
    }
    return *this;
}

CalibrationDataProvider::~CalibrationDataProvider()
{
    // Destructor implicit
}

// not used, but if the data are in a txt file, this function can be used
// will take the peaks for each source
void CalibrationDataProvider::readFromTxt(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        // std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string currentSource;
    std::vector<double> currentEnergies;

    while (std::getline(file, line))
    {
        if (line.find_first_not_of("0123456789. ") != std::string::npos)
        {
            if (!currentSource.empty())
            {
                sources.push_back(currentSource);
                energyMatrix.push_back(currentEnergies);
                currentEnergies.clear();
            }
            currentSource = line;
        }
        else
        {
            try
            {
                double energy = std::stod(line);
                currentEnergies.push_back(energy);
            }
            catch (const std::invalid_argument &)
            {
                std::cerr << "Invalid energy value: " << line << std::endl;
            }
        }
    }
    if (!currentSource.empty())
    {
        sources.push_back(currentSource);
        energyMatrix.push_back(currentEnergies);
    }

    file.close();
}

// Validates if a source exists in the loaded configuration
int CalibrationDataProvider::isSourceValid(const std::string &source)
{
    for (size_t i = 0; i < sources.size(); ++i)
    {
        if (sources[i] == source)
        {
            return i;
        }
    }
    return -1;
}

// sort the peaks in descending order
void CalibrationDataProvider::CalibrationDataProviderArray()
{
    for (auto &row : energyMatrix)
    {
        std::sort(row.begin(), row.end(), std::greater<double>());
    }
}

double *CalibrationDataProvider::getCalibratedEnergyArray(int index)
{
    if (index < 0 || index >= energyMatrix.size())
    {
        std::cerr << "Invalid index for energy array." << std::endl;
        return nullptr;
    }
    return energyMatrix[index].data();
}

int CalibrationDataProvider::getCalibratedEnergyArraySize(int index) const
{
    if (index < 0 || index >= energyMatrix.size())
    {
        std::cerr << "Invalid index for energy array." << std::endl;
        return 0;
    }
    return energyMatrix[index].size();
}

void CalibrationDataProvider::printToFile(std::ofstream &file) const
{
    for (size_t i = 0; i < sources.size(); ++i)
    {
        file << sources[i] << ": ";
        for (const auto &energy : energyMatrix[i])
        {
            file << energy << " ";
        }
        file << std::endl;
    }
}

void CalibrationDataProvider::printSources() const
{
    for (size_t i = 0; i < sources.size(); ++i)
    {
        std::cout << i << ". " << sources[i] << std::endl;
    }
}

void CalibrationDataProvider::chooseSources(int argc, char *argv[])
{
    bool dublicated = false;
    for (int i = 12; i < argc; i++)
    {
        dublicated = false;
        for (int j = 0; j < requestedSources.size(); j++)
        {
            if (requestedSources[j] == argv[i])
            {
                dublicated = true;
            }
        }
        if (!dublicated)
        {
            requestedSources.push_back(argv[i]);
        }
    }
}

void CalibrationDataProvider::chooseSources(int startPosition, int argc, char *argv[])
{
    bool dublicated = false;
    for (int i = startPosition; i < argc; i++)
    {
        dublicated = false;
        for (int j = 0; j < requestedSources.size(); j++)
        {
            if (requestedSources[j] == argv[i])
            {
                dublicated = true;
            }
        }
        if (!dublicated)
        {
            requestedSources.push_back(argv[i]);
        }
    }
}

// Gets total number of peaks across all requested sources
int CalibrationDataProvider::getNumberOfPeaks() const
{
    int totalPeaks = 0;
    for (int i = 0; i < sources.size(); i++)
    {
        for (int j = 0; j < requestedSources.size(); j++)
        {
            if (sources[i] == requestedSources[j])
            {
                totalPeaks += numberOfPeaks[i];
                break;
            }
        }
    }
    return totalPeaks;
}

int CalibrationDataProvider::getNumberOfPeaks(int position) const
{
    return numberOfPeaks[position];
}

// Creates combined energy array from requested sources
double *CalibrationDataProvider::createCalibratedSourceArray(int &size)
{
    std::string sourceNames;
    for (const auto &source : requestedSources)
    {
        sourceNames += source + " ";
    }
    ErrorHandle::getInstance().logStatus("Requested sources size: " + std::to_string(requestedSources.size()) + " Named: " + sourceNames);
    std::vector<double *> selectedEnergyArrays;
    std::vector<int> arraySizes;
    int totalSize = 0;

    for (const auto &source : requestedSources)
    {
        int index = isSourceValid(source);
        if (index != -1)
        {
            int energyArraySize = getCalibratedEnergyArraySize(index);
            if (energyArraySize > 0)
            {
                double *currentEnergyArray = getCalibratedEnergyArray(index);
                selectedEnergyArrays.push_back(currentEnergyArray);
                arraySizes.push_back(energyArraySize);
                totalSize += energyArraySize;
            }
        }
        else
        {
            std::cerr << "Invalid source name: " << source << std::endl;
            return nullptr;
        }
    }
    double *combinedEnergyArray = new double[totalSize];
    int index = 0;
    for (size_t i = 0; i < selectedEnergyArrays.size(); ++i)
    {
        int arraySize = arraySizes[i];
        std::copy(selectedEnergyArrays[i], selectedEnergyArrays[i] + arraySize, combinedEnergyArray + index);
        index += arraySize;
    }
    size = totalSize;
    return combinedEnergyArray;
}
std::string CalibrationDataProvider::cleanSourceName(const std::string &sourceName)
{
    std::string cleanedName = sourceName;
    cleanedName.erase(
        std::remove_if(cleanedName.begin(), cleanedName.end(), [](char c)
                       {
                           return c == '\"' || c == ','; // Condiția de eliminare
                       }),
        cleanedName.end());
    return cleanedName;
}

/**
 * @brief Parses JSON configuration file containing source data
 * @param filename Path to JSON file
 *
 * File format expected:
 * {
 *   "name": "sourceName",
 *   "numberOfPeaks": N,
 *   "peaks": [
 *     {"value": energy, "probability": prob}
 *   ]
 * }
 */
void CalibrationDataProvider::parseJsonFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    json j;
    file >> j;
    file.close();

    for (auto &[sourceName, sourceData] : j.items())
    {
        sources.push_back(sourceName);

        std::vector<double> energies;
        std::vector<double> probabilities;

        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Source: " << sourceName << std::endl;

        if (sourceData.contains("gammas") && sourceData["gammas"].is_object())
        {
            for (auto &[energyStr, values] : sourceData["gammas"].items())
            {
                try
                {
                    long double energyLong = std::stold(energyStr);
                    double energy = static_cast<double>(energyLong);
                    double probability = values[1]; // index 1 is the probability

                    energies.push_back(energy);
                    probabilities.push_back(probability);

                    std::cout << std::setprecision(10)
                              << "Energy: " << energy
                              << "  Probability: " << probability
                              << std::endl;
                }
                catch (const std::exception &e)
                {
                    std::cerr << "Error parsing gamma entry: " << e.what() << std::endl;
                }
            }
        }
        else
        {
            std::cerr << "Warning: 'gammas' not found or not an object in source " << sourceName << std::endl;
        }

        energyMatrix.push_back(energies);
        probabilityMatrix.push_back(probabilities);
        numberOfPeaks.push_back(static_cast<int>(energies.size()));
    }
}
