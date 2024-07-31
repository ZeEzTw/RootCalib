#include "sortEnergy.h"
#include <iostream>
#include <sstream>
#include <algorithm>

// Constructor: inițializează matricea de energie dintr-un fișier JSON
sortEnergy::sortEnergy(const std::string &filename)
{
    readFromTxt(filename);
    sortEnergyArray(); // Sortăm array-ul de energie
}

// Destructor: eliberează resursele
sortEnergy::~sortEnergy()
{
    // Destructor implicit va curăța vectorii
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
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string currentSource;
    std::vector<double> currentEnergies;

    while (std::getline(file, line))
    {
        // Verifică dacă linia este un nume de sursă (nu conține numere)
        if (line.find_first_not_of("0123456789. ") != std::string::npos)
        {
            // Dacă există o sursă curentă, o adaugă la liste
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
            // Încearcă să convertească linia într-un număr
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

    // Adaugă ultima sursă și lista de energii dacă există
    if (!currentSource.empty())
    {
        sources.push_back(currentSource);
        energyMatrix.push_back(currentEnergies);
    }

    file.close();
}
// Sortează matricea de energie în ordine descrescătoare
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
// Scrie matricea de energie sortată într-un fișier
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
