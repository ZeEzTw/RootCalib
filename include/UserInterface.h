/**
 * @brief Interactive interface for histogram calibration and source management.
 * 
 * @details
 * The `UserInterface` class provides a command-line interface for managing histograms and energy sources.
 * It is used when sources are not pre-specified via command-line arguments and influences the error messages
 * shown in the terminal while the program is running, making it perfect for debugging or testing new features.
 * If the user interface is enabled, more information will be shown in the terminal.
 * 
 * Key Features:
 * - **Interactive Source Selection**: Prompts the user to select and manage energy sources interactively.
 * - **Peak Position Modification**: Allows the user to modify peak positions in histograms after the main run is complete.
 *   This is useful for trying different calibrations, receiving data for another peak, or selecting a specific peak that was not automatically chosen.
 * - **Real-Time Feedback**: Provides immediate feedback on run status, data, calibration results, etc.
 * - **JSON Output Generation**: Generates JSON output for calibration data.
 * 
 * Methods:
 * - `askAboutPeaks`: Allows modification of peak positions in histograms.
 * - `askAboutSource`: Handles the selection and combination of energy sources.
 * - `showCalibrationInfo`: Displays calibration statistics.
 * 
 * @note This class is designed for extensibility and can be enhanced with additional user interface features as needed.
 */
#include "Histogram.h"
#include "CalibrationDataProvider.h"


class UserInterface
{
public:

    void askAboutPeaks(std::vector<Histogram> &histograms, std::ofstream &jsonFile, TFile *outputFileHistograms, TFile *outputFileCalibrated);
    double *askAboutSource(CalibrationDataProvider &energys, int &size, std::string &sourceName, int &numberOfPeaks);
    void showCalibrationInfo(const Histogram &histogram) const;

};