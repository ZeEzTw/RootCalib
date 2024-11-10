#include "../include/FileManager.h"

FileManager::FileManager(const std::string &inputFilePath, const std::string &savePath, const std::string &delila_name)
    : inputFilePath(inputFilePath), savePath(savePath), delila_name(delila_name),
      inputFile(nullptr), outputFileHistograms(nullptr),
      outputFileCalibrated(nullptr), outputFileTH2(nullptr)
{
}

FileManager::~FileManager()
{
    closeFiles();
}

void FileManager::openFiles()
{
    // Open input file
    std::cout << "inputFilePath " << inputFilePath << std::endl;
    std::cout << "savePath " << savePath << std::endl;
    inputFile = new TFile(inputFilePath.c_str(), "READ");
    std::cout << std::endl
              << "inputFile " << inputFilePath.c_str() << std::endl;
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error: Could not open input file." << std::endl;
        return;
    }

    // Get run name and base directory
    std::string runName = extractRunNumber();
    std::string baseDirectory = extractDirectoryPath();
    std::string saveDirectory;

    if (savePath.empty())
    {
        saveDirectory = baseDirectory + runName + "/";
    }
    else
    {
        saveDirectory = savePath;
        if (saveDirectory.back() != '/')
        {
            saveDirectory += '/';
        }
    }
    // Create the save directory if it doesn't exist
    if (mkdir(saveDirectory.c_str(), 0777) && errno != EEXIST)
    {
        std::cerr << "Error: Could not create save directory." << std::endl;
        return;
    }
    std::cout << "\nSave Directory: " << saveDirectory << std::endl;

    // Open JSON file
    std::string jsonFilePath = saveDirectory + runName + "_peaks_data.json";
    jsonFile.open(jsonFilePath);
    if (!jsonFile.is_open())
    {
        std::cerr << "Error: Could not open JSON file for writing." << std::endl;
        return;
    }

    // Open ROOT files for histograms
    outputFileHistograms = new TFile((saveDirectory + runName + "_peaks.root").c_str(), "RECREATE");
    outputFileCalibrated = new TFile((saveDirectory + runName + "_calibrated_histograms.root").c_str(), "RECREATE");
    outputFileTH2 = new TFile((saveDirectory + runName + "_combinedHistogram.root").c_str(), "RECREATE");

    if (!outputFileHistograms || outputFileHistograms->IsZombie())
    {
        std::cerr << "Error: Could not create output file for histograms." << std::endl;
    }
    if (!outputFileCalibrated || outputFileCalibrated->IsZombie())
    {
        std::cerr << "Error: Could not create output file for calibrated histograms." << std::endl;
    }
    if (!outputFileTH2 || outputFileTH2->IsZombie())
    {
        std::cerr << "Error: Could not create output file for combined histograms." << std::endl;
    }
}

void FileManager::closeFiles()
{
    if (inputFile)
    {
        inputFile->Close();
        delete inputFile;
        inputFile = nullptr;
    }

    if (outputFileHistograms)
    {
        outputFileHistograms->Close();
        delete outputFileHistograms;
        outputFileHistograms = nullptr;
    }

    if (outputFileCalibrated)
    {
        outputFileCalibrated->Close();
        delete outputFileCalibrated;
        outputFileCalibrated = nullptr;
    }

    if (outputFileTH2)
    {
        outputFileTH2->Close();
        delete outputFileTH2;
        outputFileTH2 = nullptr;
    }

    if (jsonFile.is_open())
    {
        jsonFile.close();
    }
}

std::string FileManager::removeFileExtension() const
{
    size_t lastDot = inputFilePath.find_last_of('.');
    if (lastDot != std::string::npos)
    {
        return inputFilePath.substr(0, lastDot);
    }
    return inputFilePath;
}

std::string FileManager::extractRunNumber() const
{
    size_t pos = 0;
    while ((pos = inputFilePath.find('_', pos)) != std::string::npos)
    {
        size_t start = pos + 1;
        size_t end = start;
        while (end < inputFilePath.size() && std::isdigit(inputFilePath[end]))
        {
            ++end;
        }
        if (end > start)
        {
            return inputFilePath.substr(start, end - start);
        }
        pos = end;
    }
    std::string fileName = removeFileExtension();
    for (char &ch : fileName)
    {
        if (ch == '/' || ch == '\\')
        {
            ch = '_';
        }
    }
    return fileName;
}

std::string FileManager::extractDirectoryPath() const
{
    size_t lastSlash = inputFilePath.find_last_of("/\\");
    if (lastSlash != std::string::npos)
    {
        return inputFilePath.substr(0, lastSlash + 1);
    }
    return "./";
}

TH2F *FileManager::getTH2Histogram() const
{
    TH2F *histogram = nullptr;
    if (inputFile)
    {
        inputFile->GetObject(delila_name.c_str(), histogram);
        if (!histogram)
        {
            std::cerr << "Error: Histogram " << delila_name << " not found in input file." << std::endl;
        }
    }
    return histogram;
}

void FileManager::saveTH2Histogram(TH2F* const th2Histogram)
{
    outputFileTH2->cd();
    th2Histogram->Write();
    // Removed redundant file closure (closed in closeFiles())
}

void FileManager::updateHistogramName(TH2F* const histogram)
{
    std::string name = histogram->GetName();
    size_t pos = name.find("_raw");
    if (pos != std::string::npos)
    {
        name.replace(pos, 4, "_calib");
        histogram->SetName(name.c_str());
    }
}
