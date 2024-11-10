#include "../include/TaskHandler.h"

int main(int argc, char *argv[])
{
    gErrorIgnoreLevel = kError; // sa scap de asta 

    ArgumentsManager argumentsManager(argc, argv);
    argumentsManager.parseJsonFile();
    argumentsManager.printArgumentsInput();
    TaskHandler taskHandler(argumentsManager);
    taskHandler.executeHistogramProcessingTask();
    std::cerr << "0" << std::endl;
    return 0;
}
