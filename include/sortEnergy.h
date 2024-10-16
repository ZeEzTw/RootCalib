#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <regex>
class sortEnergy
{
    std::vector<std::vector<double>> energyMatrix;      // Matrice de energie
    std::vector<std::vector<double>> probabilityMatrix; // Matrice de probabilitate
    std::vector<std::string> sources;                   // Numele surselor
    std::vector<int> numberOfPeaks;                     // Numărul de vâr`furi
    std::vector<std::string> requestedSources;          // Numele surselor cerute

public:
    sortEnergy(const std::string &filename);
    ~sortEnergy();

    void sortEnergyArray();
    void printToFile(std::ofstream &file) const;
    void printSources() const;
    int getSize() const { return sources.size(); }
    int getEnergyArraySize(int index) const;
    double *getEnergyArray(int index);
    void readFromTxt(const std::string &sourceLine);
    int isSourceValid(const std::string &source);
    std::string cleanSourceName(const std::string &sourceName);
    void chooseSources(int argc, char *argv[]);
    void chooseSources(int startPosition, int argc, char *argv[]);
    int getNumberOfPeaks() const;
    int getNumberOfPeaks(int position) const;
    double *createSourceArray(int &size);
    void parseJsonFile(const std::string &filename);
    const std::vector<std::string>& getRequestedSources() const { return requestedSources; }
    const std::string &getSourceName(int index) const { return sources[index]; }

private:
    // void parseJsonFile(const std::string &filename);
    // void processSource(const std::string &sourceLine, const std::string &energiesLine);
};
