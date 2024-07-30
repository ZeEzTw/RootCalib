#include "Histogram.h"
#include "Peak.h"
#include "sortEnergy.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2F.h>

// Function to open ROOT files and JSON output file
void openFiles(const char *inputFilePath, TFile *&inputFile, TFile *&outputFileHistograms, TFile *&outputFileCalibrated, std::ofstream &jsonFile) {
    inputFile = new TFile(inputFilePath, "READ");
    outputFileHistograms = new TFile("histogramsWithPeaks2.root", "RECREATE");
    outputFileCalibrated = new TFile("calibratedHistograms.root", "RECREATE");
    jsonFile.open("data.json");
}

// Function to close ROOT files and JSON output file
void closeFiles(TFile *inputFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, std::ofstream &jsonFile) {
    outputFileHistograms->Close();
    outputFileCalibrated->Close();
    inputFile->Close();
    jsonFile.close();

    delete outputFileHistograms;
    delete outputFileCalibrated;
    delete inputFile;
}

// Function to process a single 1D histogram
void processHistogram(TH1D *hist1D, int number_of_peaks, double *energyArray, int size, float Xmin, float Xmax, float FWHMmax, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated) {
    if (!hist1D || hist1D->GetMean() < 5) {
        delete hist1D;
        return;
    }

    Histogram hist(Xmin, Xmax, FWHMmax, number_of_peaks, hist1D);
    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    hist.outputPeaksDataJson(jsonFile);
    hist.printHistogramWithPeaksRoot(outputFileHistograms);
    hist.printCalibratedHistogramRoot(outputFileCalibrated);

    delete hist1D;
}

// Function to process all columns in a 2D histogram
void process2DHistogram(TH2F *h2, int number_of_peaks, double *energyArray, int size, float Xmin, float Xmax, float FWHMmax, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated) {
    if (!h2) {
        std::cerr << "Error: 2D histogram not found." << std::endl;
        return;
    }

    int number_of_columns = h2->GetNbinsX();
    for (int column = 1; column <= number_of_columns; ++column) {
        TH1D *hist1D = h2->ProjectionY(Form("hist1D_col%d", column), column, column);
        processHistogram(hist1D, number_of_peaks, energyArray, size, Xmin, Xmax, FWHMmax, jsonFile, outputFileHistograms, outputFileCalibrated);
    }
}

// Function to handle the main task of processing histograms
void processHistogramsTask(int number_of_peaks, const char *file_path, double *energyArray, int size, float Xmin, float Xmax, float FWHMmax) {
    TFile *inputFile = nullptr;
    TFile *outputFileHistograms = nullptr;
    TFile *outputFileCalibrated = nullptr;
    std::ofstream jsonFile;

    openFiles(file_path, inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);

    TH2F *h2 = nullptr;
    inputFile->GetObject("mDelila_raw", h2);
    process2DHistogram(h2, number_of_peaks, energyArray, size, Xmin, Xmax, FWHMmax, jsonFile, outputFileHistograms, outputFileCalibrated);

    closeFiles(inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);
}

// Function to sort the energy data from a file
void sortEnergyData(double *&energyArray, int &size) {
    std::ifstream inputFile("energy.txt");
    std::ofstream outputFile("output.txt");
    getEnergyArray(inputFile, energyArray, size);
    sortEnergy(energyArray, size);
    // printToFile(outputFile, energyArray, size); // Uncomment if needed
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
    sortEnergyData(energyArray, size);

    processHistogramsTask(number_of_peaks, file_path, energyArray, size, Xmin, Xmax, FWHMmax);

    delete[] energyArray;
    return 0;
}
