#include "TaskHandler.h"

TaskHandler::TaskHandler(ArgumentsManager &args)
    : argumentsManager(args), inputTH2(nullptr), energyArray(nullptr), size(0),
      fileManager(args.getHistogramFilePath(), args.getSavePath(), args.getHistogramName()) // Initialize directly in the constructor
{}

TaskHandler::~TaskHandler()
{
    if (energyArray)
    {
        delete[] energyArray;
    }
}

void TaskHandler::executeHistogramProcessingTask() {
    fileManager.openFiles();

    energyArray = initializeEnergyArray();
    if (!energyArray) {
        std::cerr << "Error: Failed to initialize energy array." << std::endl;
        return;
    }

    process2DHistogram();
    fileManager.closeFiles();
}

double* TaskHandler::initializeEnergyArray() {
    sortEnergy energyProcessor = argumentsManager.getEnergyProcessor();
    if (argumentsManager.isUserInterfaceEnabled()) {
        std::string sourcesName = argumentsManager.getSourcesName();
        int numberOfPeaks = 0;
        double* array = ui.askAboutSource(energyProcessor, size, sourcesName, numberOfPeaks);
        argumentsManager.setSourceName(sourcesName);
        argumentsManager.setNumberOfPeaks(numberOfPeaks);
        return array;
    }
    return energyProcessor.getEnergyArray(size);
}

void TaskHandler::process2DHistogram()
{
    int start_column = 0;
    inputTH2 = fileManager.getTH2Histogram(); // Use getter for inputFile
    int number_of_columns = inputTH2->GetNbinsX();

    if (argumentsManager.isDomainLimitsSet())
    {
        start_column = argumentsManager.getXminDomain();
        number_of_columns = argumentsManager.getXmaxDomain();
    }

    for (int column = start_column; column <= number_of_columns; ++column)
    {
        TH1D *hist1D = inputTH2->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D)
        {
            processSingleHistogram(hist1D);
        }
    }

    combineHistogramsIntoTH2();

    if (argumentsManager.isUserInterfaceEnabled())
    {
        ui.askAboutPeaks(histograms, fileManager.getJsonFile(), fileManager.getOutputFileHistograms(), fileManager.getOutputFileCalibrated());
    }
    // The part where UI asks if you want to change a peak
}

void TaskHandler::processSingleHistogram(TH1D* const hist1D)
{
    if (!hist1D || hist1D->GetMean() < 5)
    {
        delete hist1D;
        histograms.emplace_back();
        return;
    }

    Histogram hist;
    int histIndex = argumentsManager.getNumberColumnSpecified(histograms.size()); // Consider renaming this function
    if (argumentsManager.checkIfRunIsValid() && histIndex != -1)
    {
        hist = Histogram(
            argumentsManager.getXminFile(histIndex), argumentsManager.getXmaxFile(histIndex),
            argumentsManager.getFWHMmaxFile(histIndex), argumentsManager.getMinAmplitude(),
            argumentsManager.getMaxAmplitudeFile(histIndex), argumentsManager.getSerialFile(histIndex),
            argumentsManager.getDetTypeFile(histIndex), argumentsManager.getPolynomialFitThreshold(),
            argumentsManager.getNumberOfPeaks(), hist1D,
            argumentsManager.getHistogramNameFile(histIndex), argumentsManager.getSourcesName());
    }
    else
    {
        delete hist1D;
        histograms.emplace_back();
        return;
    }

    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    hist.outputPeaksDataJson(fileManager.getJsonFile());                      // Use getter for jsonFile
    hist.printHistogramWithPeaksRoot(fileManager.getOutputFileHistograms());  // Use getter for outputFileHistograms
    hist.printCalibratedHistogramRoot(fileManager.getOutputFileCalibrated()); // Use getter for outputFileCalibrated
    histograms.push_back(hist);

    if (argumentsManager.isUserInterfaceEnabled())
    {
        ui.showCalibrationInfo(hist);
    }

    delete hist1D;
}

void TaskHandler::combineHistogramsIntoTH2() {
    fileManager.updateHistogramName(inputTH2);
    fillTH2FromHistograms();
    fileManager.saveTH2Histogram(inputTH2);
}

void TaskHandler::fillTH2FromHistograms() {
    if (!inputTH2) return;

    inputTH2->Reset();
    int numCols = inputTH2->GetNbinsX();
    int numRows = inputTH2->GetNbinsY();

    for (size_t i = 0; i < histograms.size(); ++i) {
        TH1D *hist1D = histograms[i].getCalibratedHist();
        if (hist1D) {
            for (int binY = 1; binY <= hist1D->GetNbinsX(); ++binY) {
                double content = hist1D->GetBinContent(binY);
                if (i + 1 <= numCols && binY <= numRows) {
                    inputTH2->SetBinContent(i + 1, binY, content);
                }
            }
        }
    }
}


