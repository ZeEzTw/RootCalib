/**
 @class TaskHandler
 * @brief Coordinates and executes all necessary tasks for the main program.
 *
 * Manages file operations, user interface interactions, and histogram processing.
 * The TaskHandler class is responsible for coordinating all tasks within the program.
 * 
 * @method executeHistogramProcessingTask Executes the main task.
 * @method initializeEnergyArray Initializes the energy array with calibrated sources.
 * @method process2DHistogram Processes all histograms.
 * @method processSingleHistogram Processes a single histogram.
 * @method combineHistogramsIntoTH2 Combines histograms into a 2D histogram.
 * @method fillTH2FromHistograms Fills the 2D histogram from 1D histograms.
 */

#ifndef TASKHANDLER_H
#define TASKHANDLER_H

#include "FileManager.h"
#include "Histogram.h"
#include "Peak.h"
#include "UserInterface.h"
#include "ArgumentsManager.h"
#include <TError.h>
#include <vector>

class TaskHandler
{
private:
    ArgumentsManager &argumentsManager;
    FileManager fileManager;
    UserInterface ui;
    double *energyArray;
    int size;
    TH2F *inputTH2;
    std::vector<Histogram> histograms;

public:
    TaskHandler(ArgumentsManager &args);
    ~TaskHandler();
    void executeHistogramProcessingTask();

private:
    double *initializeEnergyArray();
    void process2DHistogram();
    void processSingleHistogram(TH1D *const hist1D);
    void combineHistogramsIntoTH2();
    void fillTH2FromHistograms();
};

#endif // TASKHANDLER_H
