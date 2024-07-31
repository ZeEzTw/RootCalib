#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
class sortEnergy
{
    std::vector<std::vector<double>> energyMatrix; // Matrice de energie
    std::vector<std::string> sources;              // Numele surselor

public:
    sortEnergy(const std::string &filename);
    ~sortEnergy();

    void sortEnergyArray();
    void printToFile(std::ofstream &file) const;
    void printSources() const;
    int getSize() const { return sources.size(); }
    int getEnergyArraySize(int index) const;
    double *getEnergyArray(int index); // Declarația corectă a metodei
    void readFromTxt(const std::string &sourceLine);

private:
    // void parseJsonFile(const std::string &filename);
    // void processSource(const std::string &sourceLine, const std::string &energiesLine);
};
