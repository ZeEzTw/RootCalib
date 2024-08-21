# Histogram Peak Extraction and Calibration Tool

This project provides tools for extracting and analyzing peaks from ROOT histograms, calibrating them using multiple sources, and saving the results. It also supports making adjustments to peaks after initial processing.

## Features

- **Peak Extraction**: Identifies the top peaks in histograms from ROOT files.
- **Filtering**: Filters peaks based on Xmin, Xmax, and FWHMmax.
- **Calibration**: Calibrates each histogram individually using specified sources.
- **File Output**: Saves the histogram with peaks, the calibrated histogram, and detailed data in a JSON file.

## Example of `data.json`

```json
{
    "Source": "Europium",
    "Histogram": "hist1D_col102",
    "NumberOfPeaks": 10,
    "Peaks": [
        {
            "Number_Peak": 1,
            "position": 91.9006,
            "amplitude": 108243,
            "sigma": 0.797756,
            "area": 191502,
            "leftLimit": 89.5073,
            "rightLimit": 94.2938
        },
       "// Additional peaks..."
    ]
}
```
Extra Features

After processing, the program offers the option to adjust peaks:

    Adjust Peak Position: You can change the position of any specified peak (e.g., move peak 5 to position 500).
    Update JSON: The data.json file is updated with the new peak positions, while preserving the old data.
    Redraw Histograms: A new spectrum is generated with the updated peak positions, while previous histograms are preserved.

This feature allows you to refine peak positions or calibrate histograms with different sets of peaks.


Installation and Usage
Clone Repository

Clone the repository using: git clone https://github.com/ZeEzTw/Eli-Europiu.git


Compilation
Compile the code with: g++ MainApp.cpp Histogram.cpp Peak.cpp sortEnergy.cpp UserInterface.cpp $(root-config --glibs --cflags --libs) -o task


Running the Program
Run the application with: ./task <number_of_peaks> <source_name> <histogram_file_path> <energy_file_path> <Xmin> <Xmax> <FWHMmax>
Example: ./task 10 Europium data.root energy.txt 10.0 1500.0 20

Argument Explanation

    <number_of_peaks>: Number of peaks to identify.
    <source_name>: Name of the source for calibration.
    <histogram_file_path>: Path to the ROOT file containing the histogram.
    <energy_file_path>: Path to the energy data file.
    <Xmin>: Minimum X limit for peak identification.
    <Xmax>: Maximum X limit for peak identification.
    <FWHMmax>: Maximum FWHM for filtering peaks.



