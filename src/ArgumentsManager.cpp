#include "../include/ArgumentsManager.h"
#include <iostream>
#include <cstring>  // pentru strcmp

ArgumentsManager::ArgumentsManager(int argc, char *argv[]) 
    : energyProcessor(energyFilePath) {  // Initialize energyProcessor with the energyFilePath
    parseArguments(argc, argv);
}

void ArgumentsManager::parseArguments(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-hf" || arg == "--histogram_file") {
            histogramFilePath = argv[++i];
        } else if (arg == "-hn" || arg == "--histogram_name") {
            TH2histogram_name = argv[++i];
        } else if (arg == "-ef" || arg == "--energy_file") {
            energyFilePath = argv[++i];
            energyProcessor = sortEnergy(energyFilePath);  // Update energyProcessor if the energy file changes
        } else if (arg == "-limits") {
            Xmin = std::stof(argv[++i]);
            Xmax = std::stof(argv[++i]);
            MinAmplitude = std::stof(argv[++i]);
            MaxAmplitude = std::stof(argv[++i]);
            FWHMmax = std::stof(argv[++i]);
        } else if (arg == "-sp" || arg == "--save_path") {
            savePath = argv[++i];
        } else if (arg == "-sources") {
            //std::cout<<"Aici"<<std::endl;
            userInterfaceStatus = false;  // Dacă se specifică sursele, UI-ul se oprește
            energyProcessor.chooseSources(i + 1, argc, argv);  // Transmit de la poziția i+1 până la final
            number_of_peaks = energyProcessor.getNumberOfPeaks();
            //std::cout<<"start: "<<i<<std::endl;
            //std::cout<<"argc: "<<argc<<std::endl;
            //std::cout<<"Sources: "<<argv[i + 1]<<std::endl;
            break;  // Nu mai continuăm, deoarece sursele sunt ultimele argumente
        }
    }
}

std::string ArgumentsManager::getSourcesName() {
    std::string sourcesName;
    const auto& usedSources = energyProcessor.getRequestedSources();  // Accesează sursele corect
    for (size_t i = 0; i < usedSources.size(); i++) {
        if (i != 0) {
            sourcesName += "+";
        }
        sourcesName += usedSources[i];
    }
    return sourcesName;
}

void ArgumentsManager::printAllArguments() {
    std::cout << "Number of peaks: " << number_of_peaks << std::endl;
    std::cout << "Histogram file path: " << histogramFilePath << std::endl;
    std::cout << "Histogram name: " << TH2histogram_name << std::endl;
    std::cout << "Energy file path: " << energyFilePath << std::endl;
    std::cout << "Xmin: " << Xmin << std::endl;
    std::cout << "Xmax: " << Xmax << std::endl;
    std::cout << "FWHMmax: " << FWHMmax << std::endl;
    std::cout << "Min amplitude: " << MinAmplitude << std::endl;
    std::cout << "Max amplitude: " << MaxAmplitude << std::endl;
    std::cout << "Save path: " << savePath << std::endl;
    std::cout << "Sources: " << getSourcesName() << std::endl;
}
