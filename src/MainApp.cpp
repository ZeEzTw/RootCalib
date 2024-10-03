#include "../include/Histogram.h"
#include "../include/Peak.h"
#include "../include/UserInterface.h"
#include "../include/ArgumentsManager.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TError.h> // Include pentru gErrorIgnoreLevel

// Function to open ROOT files and JSON output file
std::string removeFileExtension(const std::string &filePath)
{
    size_t lastDot = filePath.find_last_of('.');
    if (lastDot != std::string::npos)
    {
        return filePath.substr(0, lastDot);
    }
    return filePath;
}

void openFiles(const char *inputFilePath, TFile *&inputFile, TFile *&outputFileHistograms, TFile *&outputFileCalibrated, std::ofstream &jsonFile, TFile *&outputFileTH2, const std::string &savePath, const std::string &histogramFilePath)
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
    std::string fileName = histogramFilePath;

    // Înlăturăm extensia dacă există
    fileName = removeFileExtension(fileName);
    for (char &ch : fileName)
    {
        if (ch == '/' || ch == '\\')
        {
            ch = '_';
        }
    }
    // Construim calea completă pentru fișierul JSON
    std::string jsonFilePath = saveDirectory + fileName + ".json";
    //std::cout << "JSON file path: " << jsonFilePath << std::endl;

    jsonFile.open((jsonFilePath).c_str());
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
void convertHistogramsToTH2(const std::vector<Histogram> &histograms, TH2F *inputTH2, TFile *outputFileTH2)
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
                int xBinIndex = i + 1;

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

void processHistogram(ArgumentsManager &arguments, TH1D *hist1D, double *energyArray, int size, std::vector<Histogram> &histograms, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, UserInterface &ui)
{
    if (!hist1D || hist1D->GetMean() < 5)
    {
        delete hist1D;
        histograms.emplace_back();
        return;
    }
    Histogram hist(arguments.getXmin(), arguments.getXmax(), arguments.getFWHMmax(), arguments.getMinAmplitude(), arguments.getMaxAmplitude(), arguments.getNumberOfPeaks(), hist1D, arguments.getHistogramName(), arguments.getSourcesName());
    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    hist.outputPeaksDataJson(jsonFile);
    hist.printHistogramWithPeaksRoot(outputFileHistograms);
    hist.printCalibratedHistogramRoot(outputFileCalibrated);
    histograms.push_back(hist);

    if (arguments.isUserInterfaceEnabled())
    {
        ui.showCalibrationInfo(hist);
        std::cout << "Histograms processed: " << histograms.size() << std::endl;
    }

    delete hist1D;
}

void process2DHistogram(ArgumentsManager &arguments, TH2F *h2, double *energyArray, int size, UserInterface &ui, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated, TFile *outputFileTH2)
{
    if (!h2)
    {
        std::cerr << "4" << std::endl;
        return;
    }
    std::vector<Histogram> histograms;
    int number_of_columns = h2->GetNbinsX();
    TH2F *coppiedTh2 = (TH2F *)h2->Clone("coppiedTh2");
    for (int column = 0; column <= number_of_columns; ++column)
    {
        TH1D *hist1D = h2->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D)
        {
            processHistogram(arguments, hist1D, energyArray, size, histograms, jsonFile, outputFileHistograms, outputFileCalibrated, ui);
        }
    }

    convertHistogramsToTH2(histograms, h2, outputFileTH2);

    if (arguments.isUserInterfaceEnabled())
    {
        ui.askAboutPeaks(histograms, jsonFile, outputFileHistograms, outputFileCalibrated);
    }
}

void processHistogramsTask(ArgumentsManager &arguments)
{
    TFile *inputFile = nullptr;
    TFile *outputFileHistograms = nullptr;
    TFile *outputFileCalibrated = nullptr;
    TFile *outputFileTH2 = nullptr;
    std::ofstream jsonFile;

    openFiles(arguments.getHistogramFilePath().c_str(), inputFile, outputFileHistograms, outputFileCalibrated, jsonFile, outputFileTH2, arguments.getSavePath(), arguments.getHistogramFilePath());
    if (!inputFile)
    {
        std::cerr << "3" << std::endl;
        return;
    }

    TH2F *h2 = nullptr;
    inputFile->GetObject(arguments.getHistogramName().c_str(), h2);

    if (h2)
    {
        double *energyArray = nullptr;
        int size = 0;
        UserInterface ui;

        if (arguments.isUserInterfaceEnabled())
        {
            sortEnergy energyProcessor = arguments.getEnergyProcessor();
            energyArray = ui.askAboutSource(energyProcessor, size);
            arguments.setNumberOfPeaks(size);
            if (!energyArray)
            {
                std::cerr << "Error: Failed to retrieve energy array from UserInterface." << std::endl;
                return;
            }
        }
        else
        {
            energyArray = arguments.getEnergyProcessor().createSourceArray(size);
            if (!energyArray)
            {
                std::cerr << "2" << std::endl;
            }
        }
        process2DHistogram(arguments, h2, energyArray, size, ui, jsonFile, outputFileHistograms, outputFileCalibrated, outputFileTH2);
        delete[] energyArray;
    }
    else
    {
        std::cerr << "Error: 2D histogram '" << arguments.getHistogramName() << "' not found in the input file." << std::endl;
    }

    closeFiles(inputFile, outputFileHistograms, outputFileCalibrated, jsonFile);
}

int main(int argc, char *argv[])
{
    gErrorIgnoreLevel = kError;

    /*if (argc < 12)
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
    */
    ArgumentsManager argumentsManager(argc, argv);
    // argumentsManager.printAllArguments();
    // if (argc > 12)
    //{
    processHistogramsTask(argumentsManager);
    //}
    // else
    //{
    // processHistogramsTask(number_of_peaks, sourceName, histogramFilePath, TH2histogram_name, energyProcessor, Xmin, Xmax, FWHMmax, MinAmplitude, MaxAmplitude, savePath, histogramFilePath, true);
    //}

    std::cerr << "0" << std::endl;
    return 0;
}
