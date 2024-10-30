# Histogram Peak Extraction and Calibration Tool

This project provides tools for extracting and analyzing peaks from ROOT histograms, calibrating them using multiple sources, and saving the results. It also supports making adjustments to peaks after initial processing.

## Features

- **Peak Extraction**: Identifies the top peaks in histograms from ROOT files.
- **Filtering**: Filters peaks based on Xmin, Xmax, MinAmplitude, MaxAmplitude, and FWHMmax.
- **Calibration**: Calibrates each histogram individually using specified sources and return calibration parameters.
- **File Output**: Saves the histogram with peaks, the calibrated histogram, and detailed data in a JSON file.
- **User Interface**: Provides an optional UI for interactive calibration and adjustments.

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

    git clone https://github.com/ZeEzTw/Eli-Europiu.git

#Compilation
Compile the code with:

    g++ src/*.cpp -Iinclude $(root-config --glibs --cflags --libs) -o task


# Running the Program
 
 - 1.To run without User Interface, use this example command:

        ./task -hf "data/data.root" -hn "mDelila_raw" -ef "data/calibration_sources.json" -limits 0.0 1000000.0 0.0 1000000000.0 1000.0 -sp "output/" -detType 2 -serial "CL" -calib 1e-3 -sources 152Eu
   
   Run the program without the UI by specifying the -sources argument:

- 2.To activate User Interface, run it with the following example command:

        ./task -hf "data/data.root" -hn "mDelila_raw" -ef "data/calibration_sources.json" -limits 0.0 1000000.0 0.0 1000000000.0 1000.0 -sp "output/" -detType 2 -serial "CL" -calib 1e-3
  
- 3.Default Execution: Running the program without any additional arguments will use default values set in ArgumentsManager.h

		./task
  
Format without values puted: 
  		
    ./task -hf string -hn string -ef string -limits float float float float float -sp string -detType int -serial string -json string -calib int -domainLimits int int  -sources string

You can specify only the parameters you are interested in; any unspecified arguments will use default values or be overridden by values in the JSON file if provided. For example:

# Arguments

soruces can be puted as much as needed.


- hf / --histogram_file: Path to the histogram ROOT file. Default: data/data.root
- hn / --histogram_name: Name of the histogram within the file. Default: mDelila_raw
- ef / --energy_file: Path to the calibration source file. Default: data/calibration_sources.json
- sources: Specify calibration sources. Running without this argument starts the UI.
- limits: Limits for histogram analysis in the order: Xmin Xmax MinAmplitude MaxAmplitude FWHMmax. Default: 0.0 1000000.0 0.0 1000000000.0 1000.0
- sp / --save_path: Directory to save output files. Default: output/
- detType: Detector type used for calibration. Default: 2
- serial: Serial number of the detector. Default: CL
- json: Specify a JSON file for individual histogram configuration.
- domainLimits: Define domain limits for peak extraction in the format: xMin xMax.
- calib: Calibration polynomial fit threshold. Default: 1e-10
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



