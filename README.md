# Histogram Peak Extraction and Calibration Tool

This project provides tools for extracting and analyzing peaks from ROOT histograms, calibrating them using multiple sources, and saving the results. It also supports making adjustments to peaks after initial processing.

## Features

- **Peak Extraction**: Identifies the top peaks in histograms from ROOT files.
- **Filtering**: Filters peaks based on Xmin, Xmax, MinAmplitude, MaxAmplitude, and FWHMmax.
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
            "Number_Peak": 2,
            "position": 262.228,
            "FWHM": 1.91393,
            "area": 77169.4        
        }
    ]
    "// Additional peaks..."
}
```
## Extra Features

After processing, the program offers the option to adjust peaks:

   - Adjust Peak Position: You can change the position of any specified peak (e.g., move peak 5 to position 500).
   - Update JSON: The data.json file is updated with the new peak positions, while preserving the old data.
   - Redraw Histograms: A new spectrum is generated with the updated peak positions, while previous histograms are preserved.

This feature allows you to refine peak positions or calibrate histograms with different sets of peaks.
Is just an example to show the capability of extension with the code arhitecture.


## Installation and Usage
Clone Repository

Clone the repository using:

    git clone https://github.com/ZeEzTw/Eli-Europiu.git


# Compilation
Compile the code with:

    g++ src/*.cpp -Iinclude $(root-config --glibs --cflags --libs) -o task


# Running the Program
 - 1.To activate User Interface, run it with the following example command:

        ./task 10 Europium data/data.root mDelila_raw data/energy.txt 10.0 1500.0 20 0 10000000 output/

 - 2.To run without User Interface, use this example command:

        ./task 10 Europium data/data.root mDelila_raw data/energy.txt 10.0 1500.0 20 0 10000000 output/ Eurpoiu-152

- Run the application with: 

        ./task <number_of_peaks> <source_name> <histogram_file_path> <TH2histogram_name> <energy_file_path> <Xmin> <Xmax> <FWHMmax> <MinAmplitude> <MaxAmplitude> <save_path> <sources>

# Arguments

soruces can be puted as much as needed.

- <number_of_peaks>: Number of peaks to identify.
- <source_name>: Name of the source for calibration.
- <histogram_file_path>: Path to the ROOT file containing the histogram.
- <TH2histogram_name>: Name of the 2D histogram in the ROOT file.
- <energy_file_path>: Path to the energy data file.
- <Xmin>: Minimum X limit for peak identification.
- <Xmax>: Maximum X limit for peak identification.
- <FWHMmax>: Maximum Full Width at Half Maximum (FWHM) for filtering peaks.
- <MinAmplitude>: Minimum amplitude for peak detection.
- <MaxAmplitude>: Maximum amplitude for peak detection.
- <save_path>: Path where the results will be saved.

## Error Codes:

    0: Program finished successfully.
    1: Too few arguments to run.
    2: Source names for calibration are not valid (check spelling).
    3: Input file is not valid, or it cannot be opened (check spelling).
    4: TH2F histogram with data is not valid (check spelling).



