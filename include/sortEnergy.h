#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
class sortEnergy
{
    std::vector<std::vector<double>> energyMatrix;      // Matrice de energie
    std::vector<std::vector<double>> probabilityMatrix; // Matrice de probabilitate
    std::vector<std::string> sources;                   // Numele surselor
    std::vector<int> numberOfPeaks;                     // Numărul de vârfuri
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
    void chooseSources(int argc, char *argv[]);
    double *createSourceArray(int &size);
    void parseJsonFile(const std::string &filename);

private:
    // void parseJsonFile(const std::string &filename);
    // void processSource(const std::string &sourceLine, const std::string &energiesLine);
};
