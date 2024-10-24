#include "../include/sortEnergy.h"
#include <iostream>
#include <sstream>
#include <algorithm>

sortEnergy::sortEnergy(const std::string &filename)
{
    // readFromTxt(filename);
    parseJsonFile(filename);

    //std::ofstream outputFile("output.txt"); // Create an ofstream object
    //if (outputFile.is_open())
    //{
        //printToFile(outputFile);
        // printSources();
    //}
    // sortEnergyArray();
}

sortEnergy::~sortEnergy()
{
    // Destructor implicit
}

// Citește și parsează fișierul JSON simplificat
/*void sortEnergy::parseJsonFile(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    std::string line;
    bool inSourcesSection = false;
    std::string currentSource;
    std::vector<double> currentEnergies;
    int index = -1;

    while (std::getline(file, line)) {
        // Verifică dacă linia conține "sources"

            if (line.find("source") != std::string::npos) {
                // Salvează sursa curentă dacă există și adaugă-o în vectorul de surse

                // Extrage numele sursei
                std::string nameStart = line.find(": \"") + 3;
                std::string nameEnd = line.find("\"", nameStart);
                if (nameStart != std::string::npos && nameEnd != std::string::npos) {
                    currentSource = line.substr(nameStart, nameEnd - nameStart);
                }
                    sources.push_back(currentSource);
                    std::cout<<currentSource<<std::endl;
                    sources.push_back(currentEnergies);
                    currentEnergies.clear();
                index++;
            }
            // Verifică dacă linia conține valori de energie
            else
            {
                float energy;
                fprintf("%d\n",energy);
                energyMatrix.push_back(energy);

            }
        }
    }

// Procesează datele pentru fiecare sursă
void sortEnergy::processSource(const std::string &sourceLine, const std::string &energiesLine) {
    std::string sourceName;
    std::vector<double> energies;

    // Extrage numele sursei
    std::size_t nameStart = sourceLine.find(": \"") + 3;
    std::size_t nameEnd = sourceLine.find("\"", nameStart);
    if (nameStart != std::string::npos && nameEnd != std::string::npos) {
        sourceName = sourceLine.substr(nameStart, nameEnd - nameStart);
    }
    sources.push_back(sourceName);

    // Extrage valorile energiei
    std::stringstream ss(energiesLine);
    std::string token;
    while (std::getline(ss, token, ',')) {
        try {
            if (token.find('[') != std::string::npos) {
                token = token.substr(token.find('[') + 1);
            }
            if (token.find(']') != std::string::npos) {
                token = token.substr(0, token.find(']'));
            }
            double energy = std::stod(token);
            energies.push_back(energy);
        } catch (const std::invalid_argument&) {
            // Ignoră valorile nevalide
        }
    }

    energyMatrix.push_back(energies);
}*/
void sortEnergy::readFromTxt(const std::string &filename)
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
int sortEnergy::isSourceValid(const std::string &source)
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

void sortEnergy::sortEnergyArray()
{
    for (auto &row : energyMatrix)
    {
        std::sort(row.begin(), row.end(), std::greater<double>());
    }
}

double *sortEnergy::getEnergyArray(int index)
{
    if (index < 0 || index >= energyMatrix.size())
    {
        std::cerr << "Invalid index for energy array." << std::endl;
        return nullptr;
    }
    return energyMatrix[index].data();
}
int sortEnergy::getEnergyArraySize(int index) const
{
    if (index < 0 || index >= energyMatrix.size())
    {
        std::cerr << "Invalid index for energy array." << std::endl;
        return 0;
    }
    return energyMatrix[index].size();
}

void sortEnergy::printToFile(std::ofstream &file) const
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

void sortEnergy::printSources() const
{
    for (size_t i = 0; i < sources.size(); ++i)
    {
        std::cout << i << ". " << sources[i] << std::endl;
    }
}

void sortEnergy::chooseSources(int argc, char *argv[])
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

void sortEnergy::chooseSources(int startPosition, int argc, char *argv[])
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

int sortEnergy::getNumberOfPeaks() const
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

int sortEnergy::getNumberOfPeaks(int position) const
{
    return numberOfPeaks[position];
}
double *sortEnergy::createSourceArray(int &size)
{
    std::vector<double *> selectedEnergyArrays;
    std::vector<int> arraySizes;
    int totalSize = 0;

    for (const auto &source : requestedSources)
    {
        int index = isSourceValid(source);
        if (index != -1)
        {
            int energyArraySize = getEnergyArraySize(index);
            if (energyArraySize > 0)
            {
                double *currentEnergyArray = getEnergyArray(index);
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
std::string sortEnergy::cleanSourceName(const std::string &sourceName)
{
    std::regex pattern("[\",]");
    return std::regex_replace(sourceName, pattern, "");
}
void sortEnergy::parseJsonFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        //std::cerr << "Could not open file: " << filename << std::endl;
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
