#include "Histogram.h"
#include "Peak.h"
#include "UserInterface.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TError.h> // Include pentru gErrorIgnoreLevel 
// Function to open ROOT files and JSON output file
void openFiles(const char *inputFilePath, TFile *&inputFile, TFile *&outputFileHistograms, TFile *&outputFileCalibrated, std::ofstream &jsonFile)
{
    inputFile = new TFile(inputFilePath, "READ");
    outputFileHistograms = new TFile("histogramsWithPeaks2.root", "RECREATE");
    outputFileCalibrated = new TFile("calibratedHistograms.root", "RECREATE");
    jsonFile.open("data.json");
}

// Function to close ROOT files and JSON output file
void closeFiles(TFile *inputFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, std::ofstream &jsonFile)
{
    outputFileHistograms->Close();
    outputFileCalibrated->Close();
    inputFile->Close();
    jsonFile.close();

    delete outputFileHistograms;
    delete outputFileCalibrated;
    delete inputFile;
}

// Function to process a single 1D histogram
// Function to process a single 1D histogram
void processHistogram(TH1D *hist1D, int number_of_peaks, double *energyArray, int size, std::vector<Histogram> &histograms, float Xmin, float Xmax, float FWHMmax, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, UserInterface &ui)
{
    if (!hist1D)
    {
        std::cerr << "Error: TH1D pointer is null." << std::endl;
        return;
    }
    if (hist1D->GetMean() < 5)
    {
        delete hist1D; // Eliberare corectă a memoriei
        return;
    }
    Histogram hist(Xmin, Xmax, FWHMmax, number_of_peaks, hist1D);
    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    ui.showCalibrationInfo(hist);
    hist.outputPeaksDataJson(jsonFile);
    hist.printHistogramWithPeaksRoot(outputFileHistograms);
    hist.printCalibratedHistogramRoot(outputFileCalibrated);
    histograms.push_back(hist); // Adăugare histogramă în vectorul de histograme
    std::cout << "Histograms processed: " << histograms.size() << std::endl;
    /*
    delete hist1D; // Eliberare corectă a memoriei
    */
}

// Function to process all columns in a 2D histogram
void process2FHistogram(TH2F *h2, int number_of_peaks, double *energyArray, int size, UserInterface &ui, float Xmin, float Xmax, float FWHMmax, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated)
{
    if (!h2)
    {
        std::cerr << "Error: 2D histogram not found." << std::endl;
        return;
    }

    std::vector<Histogram> histograms;
    int number_of_columns = h2->GetNbinsX();
    std::cout << "Processing columns..." << std::endl;

    for (int column = 1; column <= number_of_columns; ++column)
    {
        TH1D *hist1D = h2->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D)
        {
            processHistogram(hist1D, number_of_peaks, energyArray, size, histograms, Xmin, Xmax, FWHMmax, jsonFile, outputFileHistograms, outputFileCalibrated, ui);
        }
        else
        {
            std::cerr << "Error: Failed to project histogram for column " << column << std::endl;
        }
    }

    // Allow the user to modify peaks in histograms after processing all columns
    ui.askAboutPeaks(histograms, jsonFile, outputFileHistograms, outputFileCalibrated);
}

// Function to handle the main task of processing histograms
void processHistogramsTask(int number_of_peaks, const std::string &file_path, sortEnergy &energyProcessor, float Xmin, float Xmax, float FWHMmax)
{
  TFile *inputFile = nullptr;
    TFile *outputFileHistograms = nullptr;
    TFile *outputFileCalibrated = nullptr;
    std::ofstream jsonFile;

    UserInterface ui;

    // Declarăm dimensiunea pentru array-ul de energie
    int size = 0;
    double *energyArray = ui.askAboutSource(energyProcessor, size);
    std::cout << size << std::endl;
    if (!energyArray)
    {
        std::cerr << "Failed to retrieve energy array from UserInterface." << std::endl;
        return;
    }

    openFiles(file_path.c_str(), inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);
 
    if (inputFile)
    {
        TH2F *h2 = nullptr;
        inputFile->GetObject("mDelila_raw", h2);
        
        if (h2)
        {
            process2FHistogram(h2, number_of_peaks, energyArray, size, ui, Xmin, Xmax, FWHMmax, jsonFile, outputFileHistograms, outputFileCalibrated);
        }
        else
        {
            std::cerr << "Error: 2D histogram not found in the input file." << std::endl;
        }

        closeFiles(inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);
    }
    else
    {
        std::cerr << "Error: Input file not opened." << std::endl;
    }

    // Eliberăm memoria alocată pentru array-ul de energie
    // delete[] energyArray;
}


// Function to sort the energy data from a file

int main(int argc, char *argv[])
{
    gErrorIgnoreLevel = kError;
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <number_of_peaks> <histogram_file_path> <energy_file_path> <Xmin> <Xmax> <FWHMmax>" << std::endl;
        return 1;
    }

    int number_of_peaks = std::atoi(argv[1]);
    const std::string histogramFilePath = argv[2];
    const std::string energyFilePath = argv[3];
    float Xmin = std::atof(argv[4]);
    float Xmax = std::atof(argv[5]);
    float FWHMmax = std::atof(argv[6]);

    // Creează un obiect sortEnergy folosind fișierul cu energiile
    sortEnergy energyProcessor(energyFilePath);

    // Folosește std::ofstream pentru a scrie într-un fișier
    std::ofstream outputFile("energiesS.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open file for writing: energiesS.txt" << std::endl;
        return 1;
    }

    energyProcessor.printToFile(outputFile);
    outputFile.close();

    // Continuă cu procesarea histogramelor
    processHistogramsTask(number_of_peaks, histogramFilePath, energyProcessor, Xmin, Xmax, FWHMmax);

    return 0;
}
