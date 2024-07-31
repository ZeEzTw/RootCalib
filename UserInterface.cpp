#include "UserInterface.h"

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
            std::cout << i << ". " << histograms[i].returnNameOfHistogram() << std::endl;
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
    std::vector<double *> selectedEnergyArrays; // Vector pentru a stoca pointerii la array-urile de energie selectate
    std::vector<int> arraySizes;                // Vector pentru a stoca dimensiunile fiecărui array de energie

    bool addMoreSources = true;
    int totalSize = 0; // Variabilă pentru a acumula dimensiunea totală

    while (addMoreSources)
    {
        // Afișăm sursele disponibile
        std::cout << "Available sources:" << std::endl;
        energys.printSources();

        // Întrebăm utilizatorul despre sursa dorită
        std::cout << "Which source do you want? (Enter the source number)" << std::endl;
        int sourceNumber;
        std::cin >> sourceNumber;

        // Verificăm validitatea sursei
        if (sourceNumber < 0 || sourceNumber >= energys.getSize())
        {
            std::cerr << "Invalid source number!" << std::endl;
            continue; // Reîncepe bucla pentru a cere o altă sursă
        }

        // Obținem pointerul la array-ul de energie și dimensiunea acestuia
        double *energyArray = energys.getEnergyArray(sourceNumber);
        int arraySize = energys.getEnergyArraySize(sourceNumber); // Obținem dimensiunea specifică a array-ului

        if (energyArray != nullptr)
        {
            selectedEnergyArrays.push_back(energyArray);
            arraySizes.push_back(arraySize);
            totalSize += arraySize; // Acumulăm dimensiunea totală
        }

        // Întrebăm utilizatorul dacă dorește să adauge o altă sursă
        std::cout << "Do you want to add more sources? (Y/N)" << std::endl;
        char answer;
        std::cin >> answer;
        addMoreSources = (answer == 'Y' || answer == 'y');
    }

    // Creăm un array unic pentru a conține toate array-urile de energie
    double *combinedEnergyArray = new double[totalSize];

    // Copiem toate valorile în `combinedEnergyArray`
    int index = 0;
    for (size_t i = 0; i < selectedEnergyArrays.size(); ++i)
    {
        int arraySize = arraySizes[i];
        std::copy(selectedEnergyArrays[i], selectedEnergyArrays[i] + arraySize, combinedEnergyArray + index);
        index += arraySize; // Avansăm indexul pentru array-ul combinat
    }

    size = totalSize; // Actualizăm dimensiunea totală

    return combinedEnergyArray; // Returnăm pointerul către array-ul combinat
}
