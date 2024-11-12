/**
 * @class FileManager
 * @brief Manages file operations, including opening, closing, and manipulating file names and paths.
 *
 * The FileManager class handles all file-related operations within the project. This includes:
 * - Opening and closing files.
 * - Editing file names and changing file paths.
 * - Providing access to specific files and histograms.
 * - Saving and updating histograms.
 *
 * The class works with various file types, including ROOT files (TFile) and JSON files (std::ofstream).
 *
 * @note The class assumes that the input files are in ROOT format and uses ROOT library classes such as TFile and TH2F.
 *
 * @param inputFilePath The path to the input file.
 * @param savePath The path where output files will be saved.
 * @param delila_name The name for TH2 histogram where the data is stored.
 */
#ifndef FILEMANAGER_H
#define FILEMANAGER_H

//#include <string>
#include <fstream>
#include <TFile.h>
#include <TH2.h>

class FileManager {
    std::string inputFilePath;
    std::string savePath;
    std::string delila_name;
    TFile* inputFile;
    TFile* outputFileHistograms;
    TFile* outputFileCalibrated;
    TFile* outputFileTH2;
    std::ofstream jsonFile;

public:
    // Constructor and Destructor
    FileManager(const std::string& inputFilePath, const std::string& savePath, const std::string& delila_name);
    ~FileManager();

    // Functions for opening and closing files
    void openFiles();
    void closeFiles();

    // Getters for private members
    TH2F* getTH2Histogram() const;
    const TFile* getInputFile() const { return inputFile; }
    std::ofstream& getJsonFile() { return jsonFile; } // Remain non-const if jsonFile needs to be modified
    TFile* getOutputFileHistograms() { return outputFileHistograms; }
    TFile* getOutputFileCalibrated() { return outputFileCalibrated; }
    const TFile* getOutputFileTH2() const { return outputFileTH2; }
    const std::string getSavePath() const { return savePath; }
    // Functions for saving and updating histograms
    void saveTH2Histogram(TH2F* const th2Histogram);
    void updateHistogramName(TH2F* const histogram);

private:
    // Funcții private pentru manipularea numelui de fișier
    std::string removeFileExtension() const;
    std::string extractRunNumber() const;
    std::string extractDirectoryPath() const;
};

#endif // FILEMANAGER_H
