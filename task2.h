#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include "TGraph.h"
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;


/*
    * Fit Gaussian function
    * @param hist - histogram
    * @param peakX - peak position
    * @param rangeWidth - range width
    * @param gaus - Gaussian function

*/
void fitGaussian(TF1 *gaus, TH1D *hist, double peakX, double rangeWidth);
/*
    * Extract the background area for the Gaussian function and return it
    * @param hist - histogram
    * @param peakX - peak position
    * @param gaussian - Gaussian function
*/
float extractBackgroundAreaForGaussian(TH1D *hist, int peakBin, TF1 *gaussian);
/*
    * Calculate the area of the Gaussian function
    * @param gaussian - Gaussian function
*/
float calculateAreaOfGaussian(TF1 *gaussian);
/*
    * Return the area of the peak after subtracting the background. 
    * @param gaussian - Gaussian function
    * @param hist - histogram
    * @param peakBin - peak position
*/
float areaPeak(TF1 *gaussian, TH1D *hist, int peakBin);
/*
    * Calculate the area of the peak, using the Gaussian function and formula resolution = sigma / mean;
    * @param gaussian - Gaussian function
*/
float calculateResolution(TF1 *gaussian);
/*
    * Eliminate the peak from the histogram, this function is used to eliminate the peak from the histogram(a temp one) to avoid double counting a peak
    * @param hist - histogram
    * @param maxBin - peak position
*/
void eliminatePeak(TH1D *hist, int maxBin);
/*
    * Find the peak in the histogram, this is the function that finds the peak in the histogram
    * @param hist - histogram
    * @param numBins - number of bins
    * @param mainHist - main histogram
*/
int findPeak(TH1D *hist, int numBins, TH1D *mainHist);
/*
    * Print the results in the file
    * @param file - file
    * @param column - column
    * @param peak - peak
    * @param peak_position - peak position
    * @param area - area
    * @param resolution - resolution
*/
void pritnInFileJson(ofstream &file, int column, int peak, int peak_position, float area, float resolution);
/*
    * Task 1 is the main function that reads the file and does the calculations and prints the results in the file
    * @param number_of_peaks - number of peaks
    * @param file_path - file path
*/
void task1(int number_of_peaks, const char *file_path);
