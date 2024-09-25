#include "../include/UserInterface.h"

void UserInterface::askAboutPeaks(std::vector<Histogram> &histograms, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated)
{
    std::cout << "Do you want to change a peak? (Y/N)" << std::endl;
    char answer;
    std::cin >> answer;

    if (answer == 'Y' || answer == 'y')
    {
        std::cout << "Available histograms:" << std::endl;
        for (int i = 0; i < histograms.size(); i++)
        {
            const char *histName = histograms[i].returnNameOfHistogram();
            if (histName && std::strlen(histName) > 0)
            {
                std::cout << i << ". " << histName << std::endl;
            }
            else
            {
                std::cout << i << ". Invalid histogram" << std::endl;
            }
        }

        while (answer == 'Y' || answer == 'y')
        {
            std::cout << "In which histogram is the peak?" << std::endl;
            int histogramNumber;
            std::cin >> histogramNumber;
 
            if (histogramNumber < 0 || histogramNumber >= histograms.size())
            {
                std::cerr << "Invalid histogram number!" << std::endl;
                continue;
            }

            std::cout << "What is the number of the peak?" << std::endl;
            int peakNumber;
            std::cin >> peakNumber;

            std::cout << "What is the new position of the peak?" << std::endl;
            double newPosition;
            std::cin >> newPosition;

            histograms[histogramNumber].changePeak(peakNumber, newPosition);
            histograms[histogramNumber].outputPeaksDataJson(jsonFile);
            histograms[histogramNumber].printHistogramWithPeaksRoot(outputFileHistograms);

            std::cout << "Do you want to change another peak? (Y/N)" << std::endl;
            std::cin >> answer;
        }
    }
}

void UserInterface::showCalibrationInfo(const Histogram &histogram)
{
    std::cout << "   matched peaks: " << histogram.getpeakMatchCount() << std::endl;
}
/*void UserInterface::askAboutCalibrationAndSources(std::vector<Histogram> &histograms, sortEnergy &energys)
{
    std::cout << "Do you want to calibrate the histograms? (Y/N)" << std::endl;
    char answer;
    std::cin >> answer;

    if (answer == 'Y' || answer == 'y')
    {
        std::cout << "Available sources:" << std::endl;
        energys.printSources();

        std::cout << "Which source do you want?" << std::endl;
        int sourceNumber;
        std::cin >> sourceNumber;

        if (sourceNumber < 0 || sourceNumber >= energys.getSize())
        {
            std::cerr << "Invalid source number!" << std::endl;
            return;
        }

        for (size_t i = 0; i < histograms.size(); ++i)
        {
            histograms[i].calibratePeaks(energys.getEnergyArray(sourceNumber), energys.getSize());
        }
    }

    // Se poate alege să recalibrezi toate histogrametele chiar și fără sursa specificată
    for (size_t i = 0; i < histograms.size(); ++i)
    {
        histograms[i].calibratePeaks(); // fără parametri
    }
}
*/

double *UserInterface::askAboutSource(sortEnergy &energys, int &size)
{
    std::vector<double *> selectedEnergyArrays;
    std::vector<int> arraySizes; 

    bool addMoreSources = true;
    int totalSize = 0;

    while (addMoreSources)
    {
        // Afișăm sursele disponibile
        std::cout << "Available sources:" << std::endl;
        energys.printSources();

        // Întrebăm utilizatorul despre sursa dorită
        std::cout << "Which source do you want? (Enter the source number)" << std::endl;
        int sourceNumber;
        std::cin >> sourceNumber;

        if (sourceNumber < 0 || sourceNumber >= energys.getSize())
        {
            std::cerr << "Invalid source number!" << std::endl;
            continue; 
        }

        double *energyArray = energys.getEnergyArray(sourceNumber);
        int arraySize = energys.getEnergyArraySize(sourceNumber); 

        if (energyArray != nullptr)
        {
            selectedEnergyArrays.push_back(energyArray);
            arraySizes.push_back(arraySize);
            totalSize += arraySize;
        }

        std::cout << "Do you want to add more sources? (Y/N)" << std::endl;
        char answer;
        std::cin >> answer;
        addMoreSources = (answer == 'Y' || answer == 'y');
    }

    double *combinedEnergyArray = new double[totalSize];

    int index = 0;
    for (size_t i = 0; i < selectedEnergyArrays.size(); ++i)
    {
        int arraySize = arraySizes[i];
        std::copy(selectedEnergyArrays[i], selectedEnergyArrays[i] + arraySize, combinedEnergyArray + index);
        index += arraySize; 
    }
    size = totalSize; 
    return combinedEnergyArray;
}
