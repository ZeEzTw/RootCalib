#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

// Clasă pentru furnizarea datelor de calibrare
class CalibrationDataProvider {
    std::vector<std::vector<double>> energyMatrix;       // Matrice de energie
    std::vector<std::vector<double>> probabilityMatrix;  // Matrice de probabilitate
    std::vector<std::string> sources;                    // Numele surselor
    std::vector<int> numberOfPeaks;                      // Numărul de vârfuri
    std::vector<std::string> requestedSources;           // Numele surselor cerute

public:
    // Constructori și Destructor
    CalibrationDataProvider(const std::string &filename);
    CalibrationDataProvider &operator=(const CalibrationDataProvider &other);
    ~CalibrationDataProvider();

    // Funcții de acces și manipulare a surselor
    int getSize() const { return sources.size(); }
    const std::vector<std::string> &getRequestedSources() const { return requestedSources; }
    const std::string &getSourceName(int index) const { return sources[index]; }

    // Funcții pentru lucrul cu datele calibrate
    double *createCalibratedSourceArray(int &size);
    int getCalibratedEnergyArraySize(int index) const;
    double *getCalibratedEnergyArray(int index);

    // Funcții pentru validarea și manipularea surselor
    void readFromTxt(const std::string &sourceLine);
    int isSourceValid(const std::string &source);
    std::string cleanSourceName(const std::string &sourceName);
    void chooseSources(int argc, char *argv[]);
    void chooseSources(int startPosition, int argc, char *argv[]);

    // Funcții pentru accesarea informațiilor despre vârfuri
    int getNumberOfPeaks() const;
    int getNumberOfPeaks(int position) const;

    // Funcții de ieșire
    void printToFile(std::ofstream &file) const;
    void printSources() const;

    // Alte funcții publice
    void CalibrationDataProviderArray();

private:
        // Metode private pentru operații interne
    void parseJsonFile(const std::string &filename);
};
