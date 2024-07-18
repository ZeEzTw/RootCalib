#include "sortEnergy.h"
void getEnergyArray(std::ifstream &file, double *&energyArray, int &size) {
    size = 0;
    double energyArrayTemp[1000];
    std::string line;
    
    while (std::getline(file, line) && size < 1000) {
        energyArrayTemp[size] = std::stod(line); // Convertim linia într-un double
        size++;
    }
    
    energyArray = new double[size];
    for (int i = 0; i < size; i++) {
        energyArray[i] = energyArrayTemp[i];
    }
}

void sortEnergy(double *&energyArray, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (energyArray[j] < energyArray[j + 1]) {
                double temp = energyArray[j];
                energyArray[j] = energyArray[j + 1];
                energyArray[j + 1] = temp;
            }
        }
    }
}

void printToFile(std::ofstream &file, double *energyArray, int size) {
    for (int i = 0; i < size; i++) {
        file << energyArray[i] << std::endl;
    }
}


