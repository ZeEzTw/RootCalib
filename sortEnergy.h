#include <fstream>
#include <fstream>

void getEnergyArray(std::ifstream &file, double *&energyArray, int &size);
void sortEnergy(double *&energyArray, int size);
void printToFile(std::ofstream &file, double *energyArray, int size);


