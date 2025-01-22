#include "TaskHandler.h"
#include "../include/ErrorHandle.h"
#include <TError.h>

TaskHandler::TaskHandler(ArgumentsManager &args)
    : argumentsManager(args), inputTH2(nullptr), energyArray(nullptr), size(0),
      fileManager(args.getHistogramFilePath(), args.getSavePath(), args.getHistogramName())
{
}

TaskHandler::~TaskHandler()
{
    if (energyArray)
    {
        delete[] energyArray;
    }
}

void TaskHandler::executeHistogramProcessingTask()
{
    ErrorHandle::getInstance().setUserInterfaceActive(argumentsManager.isUserInterfaceEnabled());
    ErrorHandle::getInstance().startProgram();
    fileManager.openFiles();
    ErrorHandle::getInstance().setPathForSave(fileManager.getSavePath());
    energyArray = initializeEnergyArray();
    if(energyArray == nullptr)
    {
        ErrorHandle::getInstance().logStatus("Energy array is null STOP the Task.");
        ErrorHandle::getInstance().saveLogFile();
        return;
    }
    if (fileManager.getTH2Histogram() == nullptr)
    {
        ErrorHandle::getInstance().logStatus("TH2F histogram is null STOP the Task.");
        ErrorHandle::getInstance().saveLogFile();
        return;
    }
    process2DHistogram();

    fileManager.closeFiles();
    ErrorHandle::getInstance().saveLogFile();
}

double *TaskHandler::initializeEnergyArray()
{
    CalibrationDataProvider energyProcessor = argumentsManager.getEnergyProcessor();
    double *array = nullptr;
    if (argumentsManager.isUserInterfaceEnabled())
    {
        std::string sourcesName;
        int numberOfPeaks = 0;
        array = ui.askAboutSource(energyProcessor, size, sourcesName, numberOfPeaks);
        argumentsManager.setSourceName(sourcesName);
        argumentsManager.setNumberOfPeaks(numberOfPeaks);
    }
    else
    {
        array = energyProcessor.createCalibratedSourceArray(size);
    }
    if (!array)
    {
        ErrorHandle::getInstance().logStatus("Failed to retrieve energy array.");
        return nullptr;
    }
    ErrorHandle::getInstance().logArrayWithCalibratedValues(array, size);
    ErrorHandle::getInstance().logStatus("Energy array initialized successfully.");
    return array;
}

void TaskHandler::process2DHistogram()
{

    int start_column = 0;
    inputTH2 = fileManager.getTH2Histogram();
    int number_of_columns = inputTH2->GetNbinsX();

    if (argumentsManager.isDomainLimitsSet())
    {
        start_column = argumentsManager.getXminDomain();
        number_of_columns = argumentsManager.getXmaxDomain();
    }
    fileManager.firstDomainJson();
    for (int column = start_column; column <= number_of_columns; ++column)
    {
        TH1D *hist1D = inputTH2->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D)
        {
            processSingleHistogram(hist1D);
        }
    }
    fileManager.lastDomainJson();
    combineHistogramsIntoTH2();

    if (argumentsManager.isUserInterfaceEnabled())
    {
        ui.askAboutPeaks(histograms, fileManager.getJsonFile(), fileManager.getOutputFileHistograms(), fileManager.getOutputFileCalibrated());
    }
    // The part where UI asks if you want to change a peak
}

void TaskHandler::processSingleHistogram(TH1D *const hist1D)
{
    if (!hist1D || hist1D->GetMean() < 5)
    {
        delete hist1D;
        histograms.emplace_back();
        // if you want to check
        // ErrorHandle::getInstance().logStatus(std::string("The mean for the histogram ") + std::to_string(histograms.size()) + " is less than 5. ");
        return;
    }

    Histogram hist;
    int histIndex = argumentsManager.getNumberColumnSpecified(histograms.size());
    if (argumentsManager.checkIfRunIsValid() && histIndex != -1)
    {
        ErrorHandle::getInstance().logStatus(std::string("Histogram: ") + std::to_string(histograms.size()) + " start to be processed.");
        ErrorHandle::getInstance().logStatus("start------------------------------------------------.");
        hist = Histogram(
            argumentsManager.getXminFile(histIndex), argumentsManager.getXmaxFile(histIndex),
            argumentsManager.getFWHMmaxFile(histIndex), argumentsManager.getMinAmplitudeFile(histIndex),
            argumentsManager.getMaxAmplitude(), argumentsManager.getSerialFile(histIndex),
            argumentsManager.getDetTypeFile(histIndex), argumentsManager.getPolynomialFitThreshold(),
            argumentsManager.getNumberOfPeaks(), hist1D,
            argumentsManager.getHistogramNameFile(histIndex), argumentsManager.getSourcesName());
    }
    else
    {
        delete hist1D;
        histograms.emplace_back();
        ErrorHandle::getInstance().logStatus(std::string("Histogram: ") + std::to_string(histograms.size()) + " is not in the Lut FIle.");
        return;
    }

    hist.findPeaks();
    hist.calibratePeaks(energyArray, size);
    hist.applyXCalibration();
    hist.outputPeaksDataJson(fileManager.getJsonFile());
    hist.printHistogramWithPeaksRoot(fileManager.getOutputFileHistograms());
    hist.printCalibratedHistogramRoot(fileManager.getOutputFileCalibrated());
    histograms.push_back(hist);

    if (argumentsManager.isUserInterfaceEnabled())
    {
        ui.showCalibrationInfo(hist);
    }

    delete hist1D;
}

void TaskHandler::combineHistogramsIntoTH2()
{
    fileManager.updateHistogramName(inputTH2);
    fillTH2FromHistograms();
    fileManager.saveTH2Histogram(inputTH2);
}

void TaskHandler::fillTH2FromHistograms()
{
    if (!inputTH2)
        return;

    inputTH2->Reset();
    int numCols = inputTH2->GetNbinsX();
    int numRows = inputTH2->GetNbinsY();

    for (size_t i = 0; i < histograms.size(); ++i)
    {
        TH1D *hist1D = histograms[i].getCalibratedHist();
        if (hist1D)
        {
            for (int binY = 1; binY <= hist1D->GetNbinsX(); ++binY)
            {
                double content = hist1D->GetBinContent(binY);
                if (i + 1 <= numCols && binY <= numRows)
                {
                    inputTH2->SetBinContent(i + 1, binY, content);
                }
            }
        }
    }
}
