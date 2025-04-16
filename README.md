# Histogram Peak Extraction and Calibration Tool

This project provides tools for extracting and analyzing peaks from ROOT histograms, calibrating them using multiple sources, and saving the results. It also supports making adjustments to peaks after initial processing. The tool is modular and extensible, offering robust error handling and an optional user interface for interactive use.

## Features

- **Peak Extraction**: Identifies the top peaks in histograms from ROOT files.
- **Filtering**: Filters peaks based on Xmin, Xmax, MinAmplitude, MaxAmplitude, and FWHMmax.
- **Calibration**: Calibrates each histogram individually using specified sources and return calibration parameters.
- **File Output**: Saves the histogram with peaks, the calibrated histogram, and detailed data in a JSON file.
- **User Interface**: Provides an optional UI for interactive calibration and adjustments.
- **Error Handling**: Logs errors and execution details, highlighting critical issues for troubleshooting.
  
## Example of data output

```json
{
  "domain": 101,
  "serial": "CL29B",
  "detType": 2,
  "PT": [
    0.1864,
    0.0001
  ],
  "pol_list": [
    1.7518,
    1.302805
  ],
  "152Eu": {
    "121.779": {
      "eff": [
        0.039331128777616554,
        0.0019667635904604525
      ],
      "res": [
        2.495,
        0.001
      ],
      "pos_ch": 92.019,
      "area": [
        3052036.973,
        2214.853
      ]
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

    g++ src/*.cpp -Iinclude $(root-config --glibs --cflags --libs) -o rootcalib


# Running the Program
 
 1. Run Without User Interface

    Use the following example commands to execute the program without the UI:
    
Basic Command:

	./rootcalib -f selected_run_163_999_eliadeS1.root -rp py_calib/LUT_RECALL.json -s 22Na

-With Path Specified:

  	/home/andrei/rootcalib/rootcalib -f selected_run_163_999_eliadeS1.root -rp /home/andrei/py_calib/LUT_RECALL.json -sc 100  -ec 140 -s 22Na

Full Constraints Example:
  
 	/home/andrei/rootcalib/rootcalib -f selected_run_163_999_eliadeS1.root -rp /home/andrei/py_calib/LUT_RECALL.json -sc 100 -ec 140 -s 22Na -sp "output/" -detType 2 -serial "CL" -calib 1e-3 -limits 0.0 1000000.0 0.0 1000000000.0 1000.0

  Run the program without the UI by specifying the -sources argument:
  
  2. Run With User Interface
  
  To activate User Interface, run it with the following example command:
  
  	./rootcalib -f selected_run_163_999_eliadeS1.root -rp /home/andrei/py_calib/LUT_RECALL.json

Format without values puted: 
  		
    ./rootcalib -f string -rp string -sc int -ec int -s string -sp string -detType int -serial string -calib float -limits float float float float float

You can specify only the parameters you are interested in; any unspecified arguments will use default values or be overridden by values in the JSON file if provided. For example:

# Arguments

soruces can be puted as much as needed.


Required Arguments

    -f / --file: Path to the ROOT histogram file.
	Format: Full path (/home/andrei/data/data.root) or shorthand (e.g., 152 for data_152_S9.root).
    -rp / --recall_path: Path to the LUT JSON file (recall profile).
	Mandatory for analysis.

Optional Arguments

	-s / --sources: List of calibration sources (e.g., 22Na, 152Eu).
	-sc / --start_channel: Starting channel for analysis.
	-ec / --end_channel: Ending channel for analysis.
	-sp / --save_path: Output directory. Default: output/.
	-detType: Detector type. Default: 2.
	-serial: Detector serial number. Default: CL.
	-calib: Calibration polynomial threshold. Default: 1e-10.
	-limits: Analysis bounds:
	Xmin Xmax MinAmplitude MaxAmplitude FWHMmax.
	Default: 0.0 1000000.0 0.0 1000000000.0 1000.0.
	-domainLimits: Peak extraction bounds: xMin xMax.
 
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
	   Fired by: ArgumentsManager::parseArguments()
	   Solution: Please provide the necessary arguments. Use -h for help.
	
	2: Source names for calibration are not valid.
	   Fired by: ArgumentsManager::validateSourceNames()
	   Solution: Check the spelling of the source names or open the file calibration_sources.json in data to see the available names. You can also run the program without -s to let the User Interface suggest sources.
	
	3: Input file is invalid or cannot be opened.
	   Fired by: Histogram::loadHistogramFile()
	   Solution: Check the spelling of the input file path.
	
	4: The TH2F histogram with data is not valid.
	   Fired by: Histogram::loadHistogram()
	   Solution: Check the spelling of the histogram name (default is mDelila_raw) or ensure the input file path is correct.
	
	5: No valid peaks for calibration (number of peaks = 0).
	   Fired by: Histogram.cpp calibratePeaksByDegree()
	   Solution: Ensure that there are valid peaks in the data.
	
	6: Output file is invalid or cannot be opened.
	   Fired by: Histogram::saveHistogram(), saveCombined(), or similar functions.
	   Solution: Check the spelling of the output file path. If this path is not specified, it will be created automatically.
	
	7: No peaks available for calibration. n = 0.
	   Fired by: Histogram::getBestDegree()
	   Solution: Probably no peaks pass the conditions (check LUT file, input data from terminal or default conditions).
	   Additional: If a histogram name was provided, it will be shown in the error message.
	
	8: LUT file not found.
	   Fired by: ArgumentsManager::parseJsonFile()
	   Solution: Check the spelling of the LUT file path.



