# Histogram Peak Extraction and Calibration Tool

This project provides tools for extracting and analyzing peaks from ROOT histograms, calibrating them using multiple sources, and saving the results. It also supports making adjustments to peaks after initial processing. The tool is modular and extensible, offering robust error handling and an optional user interface for interactive use.

## Features

- **Peak Extraction**: Identifies the top peaks in histograms from ROOT files.
- **Filtering**: Filters peaks based on Xmin, Xmax, MinAmplitude, MaxAmplitude, and FWHMmax.
- **Calibration**: Calibrates each histogram individually using specified sources and return calibration parameters.
- **File Output**: Saves the histogram with peaks, the calibrated histogram, and detailed data in a JSON file.
- **User Interface**: Provides an optional UI for interactive calibration and adjustments.
- **Error Handling**: Logs errors and execution details, highlighting critical issues for troubleshooting.
  
## Example of `data.json`

```json
{
  "domain": 102,
  "serial": "CL",
  "detType": 2,
  "PT": [0.183011, 0.000410268],
  "pol_list": [1.43531, 1.30784, -3.2675e-07],
  "152Eu": {
    "121.78": {
      "res": [1.73549e-05, 6.57896e-08],
      "pos_ch": 91.9006,
      "area": [191502, 324.311]
    }
  }
}
```

## Installation and Usage
Clone Repository

Clone the repository using:

    git clone https://github.com/ZeEzTw/RootCalib.git

#Compilation
Compile the code with:

    g++ src/*.cpp -Iinclude $(root-config --glibs --cflags --libs) -o task


# Running the Program
 
 - 1.To run without User Interface, use this example command:

        ./task -hf 152 -j "LUT_RECALL_S_20240604.json" -s "152Eu"
   
- 2.With Full Constraints
  
 		./task -hf "data/data.root" -j "LUT_RECALL_S_20240604.json" -hn "mDelila_raw" -limits 0.0 1000000.0 0.0 1000000000.0 1000.0 -sp "output/" -detType 2 -serial "CL" -calib 1e-3 -s "152Eu"

  Run the program without the UI by specifying the -sources argument:
- 3.To activate User Interface, run it with the following example command:
  
  		./task -hf 152 -j "LUT_RECALL_S_20240604.json"

Format without values puted: 
  		
    ./task -hf string -hn string -limits float float float float float -sp string -detType int -serial string -json string -calib int -domainLimits int int -sources string

You can specify only the parameters you are interested in; any unspecified arguments will use default values or be overridden by values in the JSON file if provided. For example:

# Arguments

soruces can be puted as much as needed.


Required Arguments

    -hf / --histogram_file: Path to the histogram file.
        Format: Full path (data/data.root) or shorthand (e.g., 152 for data_152_S9.root).
    -j / --json: Path to the LUT JSON file. Mandatory for analysis.

Optional Arguments

    -hn / --histogram_name: Name of the histogram. Default: mDelila_raw.
    -ef / --energy_file: Calibration source file. Default: data/calibration_sources.json. taken automatically
    -sources: List of calibration sources.
    -limits: Analysis bounds: Xmin Xmax MinAmplitude MaxAmplitude FWHMmax. Default: 0.0 1000000.0 0.0 1000000000.0 1000.0.
    -sp / --save_path: Output directory. Default: output/.
    -detType: Detector type. Default: 2.
    -serial: Detector serial number. Default: CL.
    -domainLimits: Peak extraction bounds: xMin xMax.
    -calib: Calibration polynomial threshold. Default: 1e-10.
You can specify only the parameters you need; the rest will use defaults or values from the JSON file.

## Extra Features

After processing, the program offers the option to adjust peaks:

   - Adjust Peak Position: You can change the position of any specified peak (e.g., move peak 5 to position 500).
   - Update JSON: The data.json file is updated with the new peak positions, while preserving the old data.
   - Redraw Histograms: A new spectrum is generated with the updated peak positions, while previous histograms are preserved.
   - Unmatched Peaks Detection: Detects peaks that couldn’t be matched with any source from the calibration file, which may indicate the presence of "parasite peaks" or misfits in the peak matching process.
	
This feature allows you to refine peak positions or calibrate histograms with different sets of peaks.
Is just an example to show the capability of extension with the code arhitecture.
## Error Codes:

    0: Program finished successfully.
    1: Too few arguments to run.
    2: Source names for calibration are not valid (check spelling).
    3: Input file is not valid, or it cannot be opened (check spelling).
    4: TH2F histogram with data is not valid (check spelling).
    5. No valid peaks for calibration (number of peaks = 0). Histogram.cpp calibratePeaksByDegree()



