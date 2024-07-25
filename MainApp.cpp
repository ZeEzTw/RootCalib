#include "Histogram.h"
#include "Peak.h"
#include "sortEnergy.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2F.h>

void task1(int number_of_peaks, const char *file_path, double *energyArray, int size, float Xmax, float Xmin, float FWHMmax) {
    std::ofstream file("data.json");
    TFile *fr = new TFile(file_path, "READ");
    TFile *outputFileHistograms = new TFile("histogramsWithPeaks2.root", "RECREATE");
    TFile *outputFileCalibrated = new TFile("calibratedHistograms.root", "RECREATE");

    TH2F *h1;
    fr->GetObject("mDelila_raw", h1);
    int number_of_columns = h1->GetNbinsX();
    for (int column = 1; column <= number_of_columns; column++) {
        TH1D *hist1D = h1->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D->GetMean() < 5) {
            delete hist1D;
            continue;
        }

        // Crearea unui obiect Histogram
        Histogram hist(Xmin, Xmax, FWHMmax, number_of_peaks, hist1D);
        // Apelarea funcției findPeaks()
        hist.findPeaks();
        hist.calibratePeaks(energyArray, size);
        hist.applyXCalibration();
        hist.outputPeaksDataJson(file);
        hist.printHistogramWithPeaksRoot(outputFileHistograms);
        hist.printCalibratedHistogramRoot(outputFileCalibrated);
        // Curățare resurse
        delete hist1D;
    }

    // Închiderea fișierelor
    outputFileHistograms->Close();
    outputFileCalibrated->Close();
    fr->Close();
    file.close();

    // Ștergerea obiectelor TFile
    delete outputFileHistograms;
    delete outputFileCalibrated;
    delete fr;
}
void sortTxt(double *&energyArray, int &size) {
    std::ifstream inputFile("energy.txt");
    std::ofstream outputFile("output.txt");
    getEnergyArray(inputFile, energyArray, size);
    sortEnergy(energyArray, size);
    //printToFile(outputFile, energyArray, size);
}
int main(int argc, char **argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <number_of_peaks> <file_path> <Xmin> <Xmax> <FWHMmax>" << std::endl;
        return 1;
    }
    int number_of_peaks = std::atoi(argv[1]);
    const char *file_path = argv[2];
    float Xmin = std::atof(argv[3]);
    float Xmax = std::atof(argv[4]);
    float FWHMmax = std::atof(argv[5]);

    double *energyArray = nullptr;
    int size;
    sortTxt(energyArray, size);

    task1(number_of_peaks, file_path, energyArray, size, Xmax, Xmin, FWHMmax);

    delete[] energyArray;
    return 0;
}