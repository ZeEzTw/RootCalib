#include "../include/Histogram.h"
#include "../include/Peak.h"
#include "../include/UserInterface.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TError.h> // Include pentru gErrorIgnoreLevel

// Function to open ROOT files and JSON output file
void openFiles(const char *inputFilePath, TFile *&inputFile, TFile *&outputFileHistograms, TFile *&outputFileCalibrated, std::ofstream &jsonFile, TFile *&outputFileTH2, const std::string &savePath)
{
    inputFile = new TFile(inputFilePath, "READ");
    std::string saveDirectory = savePath;
    if (!saveDirectory.empty() && saveDirectory.back() != '/')
    {
        saveDirectory += '/';
    }

    outputFileHistograms = new TFile((saveDirectory + "histogramsWithPeaks.root").c_str(), "RECREATE");
    outputFileCalibrated = new TFile((saveDirectory + "calibratedHistograms.root").c_str(), "RECREATE");
    outputFileTH2 = new TFile((saveDirectory + "combinedHistogram.root").c_str(), "RECREATE");
    jsonFile.open((saveDirectory + "data.json").c_str());
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
// Function to process all columns in a 2D histogram
void convertHistogramsToTH2(const std::vector<Histogram> &histograms, TH2F *inputTH2, TFile *outputFileTH2, int xOffset = 102)
{
    if (histograms.empty() || !inputTH2)
    {
        std::cerr << "5" << std::endl;
        return;
    }

    int number_of_columns = inputTH2->GetNbinsX();
    int number_of_binsY = inputTH2->GetNbinsY();

    inputTH2->Reset();

    int number_of_histograms = histograms.size();

    for (int i = 0; i < number_of_histograms; ++i)
    {
        TH1D *hist1D = histograms[i].getCalibratedHist();
        if (hist1D)
        {
            for (int binY = 1; binY <= hist1D->GetNbinsX(); ++binY)
            {
                double content = hist1D->GetBinContent(binY);
                int xBinIndex = i + 1 + xOffset;

                if (xBinIndex <= number_of_columns && binY <= number_of_binsY)
                {
                    inputTH2->SetBinContent(xBinIndex, binY, content);
                }
            }
        }
    }

    // Write to the TH2F object
    inputTH2->Write();

    // Write to the output file if it's valid
    if (outputFileTH2)
    {
        outputFileTH2->cd();    // Change to the output file
        inputTH2->Write();      // Write the histogram to the file
        outputFileTH2->Close(); // Close the output file
        delete outputFileTH2;   // Free the memory for the output file pointer
    }
}

void processHistogram(TH1D *hist1D, const std::string &sourceName, int number_of_peaks, double *energyArray, int size, std::vector<Histogram> &histograms, float Xmin, float Xmax, float FWHMmax, float MinAmplitude, float MaxAmplitude, const std::string &TH2histogram_name, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, UserInterface &ui, bool userInterfaceStatus)
{
    if (!hist1D || hist1D->GetMean() < 5)
    {
        delete hist1D;
        return;
    }

    Histogram hist(Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, number_of_peaks, hist1D, TH2histogram_name, sourceName);
    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    hist.outputPeaksDataJson(jsonFile);
    hist.printHistogramWithPeaksRoot(outputFileHistograms);
    hist.printCalibratedHistogramRoot(outputFileCalibrated);
    histograms.push_back(hist);

    if (userInterfaceStatus)
    {
        ui.showCalibrationInfo(hist);
        std::cout << "Histograms processed: " << histograms.size() << std::endl;
    }

    delete hist1D;
}

void process2DHistogram(TH2F *h2, const std::string &sourceName, int number_of_peaks, double *energyArray, int size, UserInterface &ui, float Xmin, float Xmax, float FWHMmax, float MinAmplitude, float MaxAmplitude, const std::string &TH2histogram_name, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, TFile *outputFileTH2, bool userInterfaceStatus)
{
    if (!h2)
    {
        std::cerr << "4" << std::endl;
        return;
    }

    std::vector<Histogram> histograms;
    int number_of_columns = h2->GetNbinsX();
    TH2F *coppiedTh2 = (TH2F *)h2->Clone("coppiedTh2");
    for (int column = 1; column <= number_of_columns; ++column)
    {
        TH1D *hist1D = h2->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D)
        {
            processHistogram(hist1D, sourceName, number_of_peaks, energyArray, size, histograms, Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, TH2histogram_name, jsonFile, outputFileHistograms, outputFileCalibrated, ui, userInterfaceStatus);
        }
    }

    convertHistogramsToTH2(histograms, h2, outputFileTH2);

    if (userInterfaceStatus)
    {
        ui.askAboutPeaks(histograms, jsonFile, outputFileHistograms, outputFileCalibrated);
    }
}

void processHistogramsTask(int number_of_peaks, const std::string &sourceName, const std::string &filePath, const std::string &TH2histogram_name, sortEnergy &energyProcessor, float Xmin, float Xmax, float FWHMmax, float MinAmplitude, float MaxAmplitude, std::string savePath, bool userInterfaceStatus)
{
    TFile *inputFile = nullptr;
    TFile *outputFileHistograms = nullptr;
    TFile *outputFileCalibrated = nullptr;
    TFile *outputFileTH2 = nullptr;
    std::ofstream jsonFile;

    openFiles(filePath.c_str(), inputFile, outputFileHistograms, outputFileCalibrated, jsonFile, outputFileTH2, savePath);

    if (!inputFile)
    {
        std::cerr << "3" << std::endl;
        return;
    }

    TH2F *h2 = nullptr;
    inputFile->GetObject(TH2histogram_name.c_str(), h2);

    if (h2)
    {
        double *energyArray = nullptr;
        int size = 0;
        UserInterface ui;

        if (userInterfaceStatus)
        {
            energyArray = ui.askAboutSource(energyProcessor, size);
            if (!energyArray)
            {
                std::cerr << "Error: Failed to retrieve energy array from UserInterface." << std::endl;
                return;
            }
        }
        else
        {
            energyArray = energyProcessor.createSourceArray(size);
            if (!energyArray)
            {
                std::cerr << "2" << std::endl;
            }
        }

        process2DHistogram(h2, sourceName, number_of_peaks, energyArray, size, ui, Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, TH2histogram_name, jsonFile, outputFileHistograms, outputFileCalibrated, outputFileTH2, userInterfaceStatus);
        delete[] energyArray;
    }
    else
    {
        std::cerr << "Error: 2D histogram '" << TH2histogram_name << "' not found in the input file." << std::endl;
    }

    closeFiles(inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);
}

int main(int argc, char *argv[])
{
    gErrorIgnoreLevel = kError;

    if (argc < 12)
    {
        std::cerr << "Usage: " << argv[0] << " <number_of_peaks> <source_name> <histogram_file_path> <TH2histogram_name> <energy_file_path> <Xmin> <Xmax> <FWHMmax> <MinAmplitude> <MaxAmplitude> <save_path> <sources> " << std::endl;
        std::cerr << "1" << std::endl;
        return 1;
    }

    int number_of_peaks = std::stoi(argv[1]);
    std::string sourceName = argv[2];
    std::string histogramFilePath = argv[3];
    std::string TH2histogram_name = argv[4];
    std::string energyFilePath = argv[5];
    float Xmin = std::stof(argv[6]);
    float Xmax = std::stof(argv[7]);
    float FWHMmax = std::stof(argv[8]);
    float MinAmplitude = std::stof(argv[9]);
    float MaxAmplitude = std::stof(argv[10]);
    std::string savePath = argv[11];

    sortEnergy energyProcessor(energyFilePath);

    if (argc > 12)
    {
        energyProcessor.chooseSources(argc, argv);
        processHistogramsTask(number_of_peaks, sourceName, histogramFilePath, TH2histogram_name, energyProcessor, Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, savePath, false);
    }
    else
    {
        processHistogramsTask(number_of_peaks, sourceName, histogramFilePath, TH2histogram_name, energyProcessor, Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, savePath, true);
    }
    std::cerr << "0" << std::endl;
    return 0;
}
