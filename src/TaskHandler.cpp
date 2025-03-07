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

//start function
void TaskHandler::executeHistogramProcessingTask()
{
    ErrorHandle::getInstance().setUserInterfaceActive(argumentsManager.isUserInterfaceEnabled());
    ErrorHandle::getInstance().startProgram();
    fileManager.openFiles();
    ErrorHandle::getInstance().setPathForSave(fileManager.getSavePath());
    energyArray = initializeEnergyArray();
    if (energyArray == nullptr)
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

//create the array with the energy values known for the specified source
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
    inputTH2 = fileManager.getTH2Histogram();
    int number_of_columns = inputTH2->GetNbinsX();
    int start_column = 0;
    int end_column = number_of_columns;

    // Set domain limits if specified
    if (argumentsManager.isDomainLimitsSet())
    {
        start_column = argumentsManager.getXminDomain();
        end_column = argumentsManager.getXmaxDomain();
        std::cout<<"Start column: "<<start_column<<std::endl;
        if (start_column < 0) start_column = 0;
        if (end_column > number_of_columns) end_column = number_of_columns;
        
        std::string logMessage = "Domain limits set - Processing from " + 
                               std::to_string(start_column) + " to " + 
                               std::to_string(end_column);
        ErrorHandle::getInstance().logStatus(logMessage);
    }

    fileManager.firstDomainJson();
    
    // First pass: Add empty histograms before start_column
    for (int column = 0; column < start_column; column++) {
        histograms.emplace_back();
    }

    // Second pass: Process histograms within range
    for (int column = start_column; column <= end_column; column++)
    {
        TH1D *hist1D = inputTH2->ProjectionY(Form("hist1D_col%d", column - 1), column, column);
        std::cout<<"hist1d_name: "<<hist1D->GetName()<<std::endl;
        if (hist1D)
        {
            std::cout<<"Processing column: "<<column<<std::endl;
            processSingleHistogram(hist1D);
        }
        else 
        {
            histograms.emplace_back();
        }
    }

    // Third pass: Add empty histograms after end_column
    for (int column = end_column + 1; column <= number_of_columns; ++column) {
        histograms.emplace_back();
    }

    fileManager.lastDomainJson();
    combineHistogramsIntoTH2();
    
    if (argumentsManager.isUserInterfaceEnabled())
    {
        ui.askAboutPeaks(histograms, fileManager.getJsonFile(), 
                        fileManager.getOutputFileHistograms(), 
                        fileManager.getOutputFileCalibrated());
    }
}


void TaskHandler::processSingleHistogram(TH1D *const hist1D)
{
    if (!hist1D || hist1D->GetMean() < 5)
    {
        if(hist1D->GetMean() < 5)
        std::cout<<"Skipped1"<<std::endl;
        else
        std::cout<<"Skipped2"<<std::endl;
        delete hist1D;
        histograms.emplace_back();
        // if you want to check
        // ErrorHandle::getInstance().logStatus(std::string("The mean for the histogram ") + std::to_string(histograms.size()) + " is less than 5. ");
        return;
    }

    Histogram hist;
    int histIndex = argumentsManager.getNumberColumnSpecified(histograms.size());
    ErrorHandle::getInstance().logStatus(std::string("Histogram: ") + std::to_string(histograms.size()) + " start to be processed.");
    ErrorHandle::getInstance().logStatus("start------------------------------------------------.");
    hist = Histogram(
        argumentsManager.getXminFile(histIndex), argumentsManager.getXmaxFile(histIndex),
        argumentsManager.getFWHMmaxFile(histIndex), argumentsManager.getMinAmplitudeFile(histIndex),
        argumentsManager.getMaxAmplitude(), argumentsManager.getSerialFile(histIndex),
        argumentsManager.getDetTypeFile(histIndex), argumentsManager.getPolynomialFitThreshold(),
        argumentsManager.getNumberOfPeaks(), hist1D,
        argumentsManager.getHistogramNameFile(histIndex), argumentsManager.getSourcesName());

    //all the function that take care of the task 
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
