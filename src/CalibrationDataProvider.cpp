#include "../include/CalibrationDataProvider.h"

CalibrationDataProvider::CalibrationDataProvider(const std::string &filename)
{
    parseJsonFile(filename);
}

CalibrationDataProvider& CalibrationDataProvider::operator=(const CalibrationDataProvider &other)
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

void CalibrationDataProvider::readFromTxt(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        //std::cerr << "Could not open file: " << filename << std::endl;
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
int CalibrationDataProvider::isSourceValid(const std::string &source)
{
    for (size_t i = 0; i < sources.size(); ++i)
    {
        //std::cout << sources[i] << std::endl;
        if (sources[i] == source)
        {
            return i;
        }
    }
    return -1;
}

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
double *CalibrationDataProvider::createCalibratedSourceArray(int &size)
{
    std::cout<<"requestedSources.size()"<<requestedSources.size()<<std::endl;
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
std::string CalibrationDataProvider::cleanSourceName(const std::string &sourceName) {
    std::string cleanedName = sourceName;
    cleanedName.erase(
        std::remove_if(cleanedName.begin(), cleanedName.end(), [](char c) {
            return c == '\"' || c == ',';  // Condiția de eliminare
        }),
        cleanedName.end()
    );
    return cleanedName;
}
void CalibrationDataProvider::parseJsonFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string currentSource;
    int peakCount = 0;

    while (std::getline(file, line))
    {
        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

        if (line.find("name") != std::string::npos)
        {
            currentSource = line.substr(line.find(":") + 1);
            currentSource = cleanSourceName(currentSource);
            sources.push_back(currentSource);
        }

        if (line.find("numberOfPeaks") != std::string::npos)
        {
            peakCount = std::stoi(line.substr(line.find(":") + 1));
            numberOfPeaks.push_back(peakCount);
        }

        if (line.find("peaks") != std::string::npos || line.find("Peaks") != std::string::npos)
        {
            std::vector<double> energies;
            std::vector<double> probabilities;

            while (std::getline(file, line) && line.find("]") == std::string::npos)
            {
                if (line.find("value") != std::string::npos || line.find("Energy_keV") != std::string::npos)
                {
                    double energy = std::stod(line.substr(line.find(":") + 1));
                    energies.push_back(energy);
                }
                else
                {
                    continue;
                }
                std::getline(file, line);
                if (line.find("probability") != std::string::npos || line.find("Probability") != std::string::npos)
                {
                    double probability = std::stod(line.substr(line.find(":") + 1));
                    probabilities.push_back(probability);
                }
                else
                {
                    continue;
                }
            }

            energyMatrix.push_back(energies);
            probabilityMatrix.push_back(probabilities);
        }
    }

    file.close();
}
