#include "../include/Histogram.h"
#include "../include/EliadeMathFunctions.h"
#include "../include/ErrorHandle.h"
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <limits>
//#include <TCanvas.h>
//#include <TLatex.h>
//#include <TGraph.h>

// Global constants
namespace
{
    constexpr float MAX_DISTANCE = 10;
    constexpr float MIN_DISTANCE = 1.9f;
}

// Constructor implementations
Histogram::Histogram() : xMin(0), xMax(0), maxFWHM(0), minAmplitude(0), maxAmplitude(0),
                         numberOfPeaks(0), mainHist(nullptr), tempHist(nullptr), calibratedHist(nullptr),
                         m(0), b(0), polynomialDegree(0), peakMatchCount(0), // Changed from polinomDegree
                         polynomialFitThreshold(1e-3),
                         TH2histogram_name("An empty histogram"),
                         sourceName("Empty source"),
                         peakCount(0),
                         totalArea(0),
                         totalAreaError(0)
{
}

Histogram::Histogram(int xMin, int xMax, int maxFWHM, float minAmplitude, float maxAmplitude,
                     const std::string &serial, int detType, float polynomialFitThreshold, int numberOfPeaks,
                     TH1D *mainHist, const std::string &TH2histogram_name, const std::string &sourceName)
    : xMin(xMin), xMax(xMax), maxFWHM(maxFWHM), minAmplitude(minAmplitude), maxAmplitude(maxAmplitude),
      serial(serial), detType(detType), polynomialFitThreshold(polynomialFitThreshold), numberOfPeaks(numberOfPeaks),
      mainHist(mainHist), tempHist(nullptr), calibratedHist(nullptr),
      m(0), b(0), polynomialDegree(0), peakMatchCount(0),
      TH2histogram_name(TH2histogram_name), sourceName(sourceName),
      peakCount(0), totalArea(0), totalAreaError(0)
{
    if (mainHist)
    {
        this->tempHist = (TH1D *)mainHist->Clone();
        this->calibratedHist = (TH1D *)mainHist->Clone();
    }
}

Histogram::~Histogram()
{
    // Consider using smart pointers to manage memory automatically
    // if (tempHist) { delete tempHist; tempHist = nullptr; }
    // if (calibratedHist) { delete calibratedHist; calibratedHist = nullptr; }
}

Histogram::Histogram(const Histogram &histogram)
    : xMin(histogram.xMin), xMax(histogram.xMax), maxFWHM(histogram.maxFWHM),
      minAmplitude(histogram.minAmplitude), maxAmplitude(histogram.maxAmplitude),
      numberOfPeaks(histogram.numberOfPeaks), serial(histogram.serial), detType(histogram.detType),
      polynomialFitThreshold(histogram.polynomialFitThreshold), m(histogram.m), b(histogram.b),
      peakMatchCount(histogram.peakMatchCount), polynomialDegree(histogram.polynomialDegree), // Changed from polinomDegree
      TH2histogram_name(histogram.TH2histogram_name), sourceName(histogram.sourceName),
      peakCount(histogram.peakCount), totalArea(histogram.totalArea), totalAreaError(histogram.totalAreaError)
{
    mainHist = (histogram.mainHist) ? (TH1D *)histogram.mainHist->Clone() : nullptr;
    tempHist = (histogram.tempHist) ? (TH1D *)histogram.tempHist->Clone() : nullptr;
    calibratedHist = (histogram.calibratedHist) ? (TH1D *)histogram.calibratedHist->Clone() : nullptr;
    peaks = histogram.peaks;
    coefficients = histogram.coefficients;
}

Histogram &Histogram::operator=(const Histogram &histogram)
{
    if (this != &histogram)
    {
        if (mainHist)
            delete mainHist;
        if (tempHist)
            delete tempHist;
        if (calibratedHist)
            delete calibratedHist;

        xMin = histogram.xMin;
        xMax = histogram.xMax;
        maxFWHM = histogram.maxFWHM;
        minAmplitude = histogram.minAmplitude;
        maxAmplitude = histogram.maxAmplitude;
        numberOfPeaks = histogram.numberOfPeaks;
        serial = histogram.serial;
        detType = histogram.detType;
        polynomialFitThreshold = histogram.polynomialFitThreshold;
        m = histogram.m;
        b = histogram.b;
        peakMatchCount = histogram.peakMatchCount;
        polynomialDegree = histogram.polynomialDegree;
        TH2histogram_name = histogram.TH2histogram_name;
        sourceName = histogram.sourceName;
        peakCount = histogram.peakCount;
        totalArea = histogram.totalArea;
        totalAreaError = histogram.totalAreaError;

        mainHist = (histogram.mainHist) ? (TH1D *)histogram.mainHist->Clone() : nullptr;
        tempHist = (histogram.tempHist) ? (TH1D *)histogram.tempHist->Clone() : nullptr;
        calibratedHist = (histogram.calibratedHist) ? (TH1D *)histogram.calibratedHist->Clone() : nullptr;
        peaks = histogram.peaks;
        coefficients = histogram.coefficients;
    }
    return *this;
}

// peak detection section
void Histogram::findPeaks()
{
    int result = 0;
    int count = 0;
    while (result == 0 && count < numberOfPeaks)
    {
        result = detectAndFitPeaks();
        if (result == -1)
        {
            break;
        }
        count++;
    }
    ErrorHandle::getInstance().logStatus("Peaks detected: " + std::to_string(peaks.size()));
}
void Histogram::eliminatePeak(const Peak &peak)
{
    double mean = peak.getMean();
    double sigma = peak.getSigma();
    int leftLimit = static_cast<int>(mean - 3 * sigma);
    int rightLimit = static_cast<int>(mean + 3 * sigma);

    if (leftLimit < 1)
        leftLimit = peak.getPosition() - 5;
    if (rightLimit < 1)
        rightLimit = peak.getPosition() + 5;
    if (rightLimit > tempHist->GetNbinsX())
        rightLimit = tempHist->GetNbinsX();
    for (int i = leftLimit; i <= rightLimit; ++i)
    {
        tempHist->SetBinContent(i, 0);
    }
}

bool Histogram::checkConditions(const Peak &peak) const
{
    double FWHM = peak.getFWHM();
    double peakPosition = peak.getPosition();
    bool condition1 = peakPosition <= xMax && peakPosition >= xMin;
    bool condition2 = FWHM <= maxFWHM;
    bool condition3 = peak.getAmplitude() > minAmplitude && peak.getAmplitude() < maxAmplitude;

    return condition1 && condition2 && condition3;
}

TF1 *Histogram::createGaussianFit(int maxBin)
{
    float maxPeakX = mainHist->GetXaxis()->GetBinCenter(maxBin);
    // is good to use [0]*exp(-0.5*((x-[1])/[2])**2) + [3] + ([4]*x) to fit the gaussian
    TF1 *gaus = new TF1(Form("gausFit_%d", peakCount), "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + ([4]*x) + ([5]*x*x)", maxPeakX - 10, maxPeakX + 10);
    gaus->SetParameters(tempHist->GetBinContent(maxBin), maxPeakX, 1.0, 0.0, 0.0, 0.0);
    gaus->SetParLimits(2, 0.1, 10.0);
    return gaus;
}

int Histogram::findMaxBin()
{
    double threshold = 0.001;
    float maxPeakY = 0;
    int maxBin = 0;
    double peakWithoutBackground = 0;
    double leftLimit = MIN_DISTANCE;
    double rightLimit = MIN_DISTANCE;

    for (int bin = 1; bin <= tempHist->GetNbinsX(); ++bin)
    {
        float binContent = tempHist->GetBinContent(bin);
        if (binContent == 0)
            continue;
        double leftLimityContent = mainHist->GetBinContent(mainHist->FindBin(bin - leftLimit));
        double rightLimityContent = mainHist->GetBinContent(mainHist->FindBin(bin + rightLimit));
        double peakWithoutBackgroundTemp = binContent - ((leftLimityContent + rightLimityContent) / 2);
        if (leftLimityContent < threshold)
        {
            peakWithoutBackgroundTemp = binContent - rightLimityContent;
        }
        else if (rightLimityContent < threshold)
        {
            peakWithoutBackgroundTemp = binContent - leftLimityContent;
        }
        else
        {
            peakWithoutBackgroundTemp = binContent - ((leftLimityContent + rightLimityContent) / 2);
        }
        // double peakWithoutBackgroundTemp1 = binContent - (mainHist->GetBinContent(mainHist->FindBin(bin - leftLimit)) + mainHist->GetBinContent(mainHist->FindBin(bin + rightLimit))) / 2;
        if (peakWithoutBackgroundTemp > peakWithoutBackground)
        {
            peakWithoutBackground = peakWithoutBackgroundTemp;
            maxPeakY = binContent;
            maxBin = bin;
        }
    }

    return maxBin;
}

// Calibration section
int Histogram::detectAndFitPeaks()
{
    int maxBin = findMaxBin();
    if (maxBin == 0)
    {
        return -1;
    }

    TF1 *gaus = createGaussianFit(maxBin);
    tempHist->Fit(gaus, "RQ");
    peaks.emplace_back(gaus, mainHist);

    if (!checkConditions(peaks.back()))
    {
        ErrorHandle::getInstance().logStatus("Peak " + std::to_string(peaks.back().getPosition()) + " does not meet the conditions, peak procces stops here.");
        eliminatePeak(peaks.back());
        peaks.pop_back();
        delete gaus;
        return 0;
    }

    // Validate peak quality
    if (!isValidPeak(peaks.back()))
    {
        peaks.pop_back();
        delete gaus;
        return -1;
    }

    peakCount++;
    eliminatePeak(peaks.back());

    delete gaus;
    return 0;
}

bool Histogram::isValidPeak(const Peak &peak) const
{
    return peak.getArea() > 0 &&
           (peaks.empty() || peak.getPosition() != peaks[peaks.size() - 2].getPosition()) &&
           peak.getPosition() >= 0;
}

bool Histogram::checkPredictedEnergies(double predictedEnergy, const double knownEnergies[], int size, float errorAdmitted, double &valueAssociatedWith) const
{
    double minError = std::numeric_limits<double>::max();
    for (int i = 0; i < size; ++i)
    {
        double error = std::abs(predictedEnergy - knownEnergies[i]);
        if (error < minError)
        {
            minError = error;
            valueAssociatedWith = knownEnergies[i];
        }
    }

    return minError < errorAdmitted;
}

void Histogram::calibratePeaks(const double knownEnergies[], int size)
{
    double bestM = 0.0;
    double bestB = 0.0;
    int bestCorrelation = 0;
    double valueAssociatedWith = 0.0;

    for (double m = 0.01; m <= 5.0; m += 0.0001)
    {
        std::vector<double> associatedValues(peaks.size(), 0.0);
        int correlations = 0;
        int peakCount = 0;

        for (const auto &peak : peaks)
        {
            double predictedEnergy = m * peak.getPosition() + b;

            if (checkPredictedEnergies(predictedEnergy, knownEnergies, size, 10, valueAssociatedWith))
            {
                ++correlations;
                associatedValues[peakCount] = valueAssociatedWith;
            }
            peakCount++;
        }
        if (correlations > bestCorrelation)
        {
            bestCorrelation = correlations;
            bestM = m;
            for (int i = 0; i < peaks.size(); ++i)
            {
                peaks[i].setAssociatedPosition(associatedValues[i]);
            }
        }
    }
    ErrorHandle::getInstance().logStatus("Peaks asociated with calibrated ones: " + std::to_string(bestCorrelation));
    peakMatchCount = bestCorrelation;
    calibratePeaksByDegree();
}

// getting polynomial degree + values

void Histogram::calibratePeaksByDegree()
{
    calibrationDegree = 1;

    std::vector<double> positions;
    std::vector<double> energies;

    for (const auto &peak : peaks)
    {
        if (peak.getAssociatedPosition() > 0)
        {
            positions.push_back(peak.getPosition());
            energies.push_back(peak.getAssociatedPosition());
        }
    }

    int n = positions.size();
    if (n == 0)
    {
        ErrorHandle::getInstance().errorHandle(ErrorHandle::NO_PEAKS_FOR_CALIBRATION);
        return;
    }

    // Gradul maxim bazat pe numărul de puncte (n - 1)
    int maxDegree = std::min(n - 1, 3);
    for (int currentDegree = 1; currentDegree <= maxDegree; ++currentDegree)
    {
        std::vector<std::vector<double>> X(n, std::vector<double>(currentDegree + 1));
        std::vector<double> Y(n);

        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= currentDegree; ++j)
            {
                X[i][j] = std::pow(positions[i], j);
            }
            Y[i] = energies[i];
        }

        std::vector<std::vector<double>> XtX = EliadeMathFunctions::multiplyTransposeMatrix(X);
        std::vector<double> XtY = EliadeMathFunctions::multiplyTransposeVector(X, Y);
        std::vector<double> coeffs = EliadeMathFunctions::solveSystem(XtX, XtY);

        if (std::abs(coeffs[currentDegree]) >= polynomialFitThreshold)
        {
            calibrationDegree = currentDegree;
            coefficients = coeffs;
        }
    }
}

// V2 calibration section //removed, available in the previous version on github

// Apply calibration section
void Histogram::initializeCalibratedHist()
{
    if (mainHist == nullptr)
    {
        return;
    }

    if (calibratedHist == nullptr)
    {
        int numBins = mainHist->GetNbinsX();
        double xMin = mainHist->GetXaxis()->GetXmin();
        double xMax = mainHist->GetXaxis()->GetXmax();
        calibratedHist = new TH1D("calibratedHist", "Calibrated Histogram", numBins, xMin, xMax);
    }
}

void Histogram::interpolateBins(int start_bin, int end_bin, double start_value, double end_value, double start_position, double end_position)
{
    if (end_bin > start_bin)
    {
        double distance = end_position - start_position;

        for (int i = start_bin + 1; i < end_bin; ++i)
        {
            double fraction = (calibratedHist->GetXaxis()->GetBinCenter(i) - start_position) / distance;

            double interpolated_content = start_value + fraction * (end_value - start_value);

            calibratedHist->SetBinContent(i, interpolated_content);
        }
    }
}

void Histogram::setCalibratedBinContent(int original_bin)
{
    double binCenter_original = mainHist->GetXaxis()->GetBinCenter(original_bin);

    double binCenter_calibrated = evaluateCalibrationPolynomial(binCenter_original);

    double binContent = mainHist->GetBinContent(original_bin);

    int corresponding_bin_calibrated = calibratedHist->GetXaxis()->FindBin(binCenter_calibrated);

    if (corresponding_bin_calibrated < 1 || corresponding_bin_calibrated > calibratedHist->GetNbinsX())
    {
        calibratedHist->SetBinContent(corresponding_bin_calibrated, 0);
    }
    else
    {
        calibratedHist->SetBinContent(corresponding_bin_calibrated, binContent);
    }
}

double Histogram::evaluateCalibrationPolynomial(double binCenter_original) const
{
    double binCenter_calibrated = 0;
    double binCenter_power = 1;
    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        binCenter_calibrated += coefficients[i] * binCenter_power;
        binCenter_power *= binCenter_original;
    }

    return binCenter_calibrated;
}

void Histogram::applyXCalibration()
{
    initializeCalibratedHist();
    int numBinsOriginal = mainHist->GetNbinsX();

    for (int new_bin = 1; new_bin <= numBinsOriginal; ++new_bin)
    {
        setCalibratedBinContent(new_bin);

        if (new_bin < numBinsOriginal)
        {
            double nextBinCenter_original = mainHist->GetXaxis()->GetBinCenter(new_bin + 1);
            double nextBinCenter_calibrated = evaluateCalibrationPolynomial(nextBinCenter_original);
            double nextBinContent = mainHist->GetBinContent(new_bin + 1);

            double binCenter_original_calibrated = evaluateCalibrationPolynomial(mainHist->GetXaxis()->GetBinCenter(new_bin));
            int current_bin_calibrated = calibratedHist->GetXaxis()->FindBin(binCenter_original_calibrated);
            int next_bin_calibrated = calibratedHist->GetXaxis()->FindBin(nextBinCenter_calibrated);

            interpolateBins(current_bin_calibrated, next_bin_calibrated,
                            mainHist->GetBinContent(new_bin), nextBinContent,
                            binCenter_original_calibrated, nextBinCenter_calibrated);
        }
    }
}

// output section
void Histogram::outputPeaksDataJson(std::ofstream &jsonFile)
{
    // jsonFile << "[\n";
    jsonFile << "\t{\n";
    jsonFile << "\t\t\"domain\": " << getMainHistName() << ",\n";
    jsonFile << "\t\t\"serial\": \"" << serial << "\",\n";
    jsonFile << "\t\t\"detType\": " << detType << ",\n";
    jsonFile << "\t\t\"PT\": [";
    jsonFile << getPT() << ", " << getPTError();
    jsonFile << "],\n";

    jsonFile << "\t\t\"pol_list\": [\n";
    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        jsonFile << "\t\t\t" << coefficients[i];
        if (i < coefficients.size() - 1)
        {
            jsonFile << ",";
        }
        jsonFile << "\n";
    }
    jsonFile << "\t\t],\n";

    jsonFile << "\t\t\"" << sourceName << "\": {\n";

    for (size_t i = 0; i < peaks.size(); ++i)
    {
        jsonFile << "\t\t\t\"" << peaks[i].getAssociatedPosition() << "\": {\n";
        // jsonFile << "\t\t\t\t\"eff\": [" << peaks[i].getEfficiency() << ", " << peaks[i].getEffError() << "],\n";
        jsonFile << "\t\t\t\t\"res\": [" << peaks[i].calculateResolution() << ", " << peaks[i].calculateResolutionError() << "],\n";
        jsonFile << "\t\t\t\t\"pos_ch\": " << peaks[i].getPosition() << ",\n";
        jsonFile << "\t\t\t\t\"area\": [" << peaks[i].getArea() << ", " << peaks[i].getAreaError() << "]\n";

        if (i < peaks.size() - 1)
        {
            jsonFile << "\t\t\t},\n";
        }
        else
        {
            jsonFile << "\t\t\t\n";
        }
    }

    jsonFile << "\t\t}\n";
    jsonFile << "\t},\n";
    // jsonFile << "]\n";
    ErrorHandle::getInstance().logStatus("Peaks data saved successfully in Json.");
    ErrorHandle::getInstance().logStatus("end------------------------------------------------.");

}

void Histogram::printHistogramWithPeaksRoot(TFile *outputFile)
{
    if (!outputFile || outputFile->IsZombie())
    {
        ErrorHandle::getInstance().logStatus("Error: Could not open file for writing in printHistogramWithPeaksRoot, filed output its not send corectly.");
        return;
    }

    outputFile->cd();

    mainHist->GetListOfFunctions()->Clear();

    for (int i = 0; i < peaks.size(); ++i)
    {
        double newPosition = peaks[i].getPosition();
        TF1 *gaussianFunction = createGaussianFit(newPosition);
        if (gaussianFunction)
        {
            if (peaks[i].getAssociatedPosition() == 0)
            {
                gaussianFunction->SetLineColor(kGreen);
            }
            gaussianFunction->SetName(Form("gaussian_%d", i));
            gaussianFunction->SetTitle(Form("Gaussian %d", i));
            mainHist->Fit(gaussianFunction, "RQ+");
        }
    }

    mainHist->Write();
}

void Histogram::printCalibratedHistogramRoot(TFile *outputFile) const
{
    if (!outputFile || outputFile->IsZombie())
    {
        ErrorHandle::getInstance().logStatus("Error: Could not open file for writing in printCalibratedHistogramRoot, filed output its not send corectly.");
        return;
    }
    outputFile->cd();
    calibratedHist->Write();
}

// set/get functions section
void Histogram::setTotalArea()
{
    totalArea = 0;
    for (int bin = 1; bin <= mainHist->GetNbinsX(); ++bin)
    {
        totalArea += mainHist->GetBinContent(bin) * mainHist->GetBinWidth(bin);
    }
}

void Histogram::setTotalAreaError()
{
    totalAreaError = 0;
    for (int bin = 1; bin <= mainHist->GetNbinsX(); ++bin)
    {
        double binError = mainHist->GetBinError(bin);
        double binWidth = mainHist->GetBinWidth(bin);
        totalAreaError += std::pow(binError * binWidth, 2);
    }
    totalAreaError = std::sqrt(totalAreaError);
}

float Histogram::getPT()
{
    setTotalArea();
    setTotalAreaError();
    float areaPeak = 0;
    for (const auto &peak : peaks)
    {
        areaPeak += peak.getArea();
    }
    return totalArea > 0 ? areaPeak / totalArea : 0.0f;
}

float Histogram::getPTError()
{
    float areaPeak = 0;
    float areaPeakError = 0;

    for (const auto &peak : peaks)
    {
        areaPeak += peak.getArea();
        areaPeakError += std::pow(peak.getAreaError(), 2);
    }
    areaPeakError = std::sqrt(areaPeakError);
    float pt = totalArea > 0 ? areaPeak / totalArea : 0.0f;
    float ptError = pt * std::sqrt(std::pow(areaPeakError / areaPeak, 2) + std::pow(totalAreaError / totalArea, 2));
    return ptError > 0 ? ptError : 0.0f;
}

void Histogram::findStartOfPeak(Peak &peak)
{
    double mean = peak.getMean();
    double sigma = peak.getSigma();

    double leftLimitPosition = mean - 2 * sigma;
    double rightLimitPosition = mean + 2 * sigma;

    double left = std::fabs(leftLimitPosition - peak.getPosition());
    double right = std::fabs(rightLimitPosition - peak.getPosition());

    if (left > MAX_DISTANCE || left < MIN_DISTANCE)
    {
        leftLimitPosition = peak.getPosition() - MIN_DISTANCE;
    }
    if (right > MAX_DISTANCE || right < MIN_DISTANCE)
    {
        rightLimitPosition = peak.getPosition() + MIN_DISTANCE;
    }
}

const char *Histogram::returnNameOfHistogram() const
{
    if (mainHist)
    {
        return mainHist->GetName();
    }
    else
    {
        return "Invalid histogram";
    }
}

// extra function section
void Histogram::changePeak(int peakNumber, double newPosition)
{
    if (peakNumber < 0 || peakNumber >= static_cast<int>(peaks.size()))
    {
        std::cerr << "Error: Invalid peak number." << std::endl;
        return;
    }

    // Use C++-style cast
    TF1 *gaus = createGaussianFit(static_cast<int>(newPosition));
    if (gaus == nullptr)
    {
        std::cerr << "Error: Failed to create Gaussian fit." << std::endl;
        return;
    }

    tempHist->Fit(gaus, "RQ");

    eliminatePeak(peaks[peakNumber]);

    peaks[peakNumber] = Peak(gaus, mainHist);
    std::cout << "Peak number: " << peakNumber << std::endl;
    std::cout << "Peak position: " << peaks[peakNumber].getPosition() << std::endl;
    if (!checkConditions(peaks[peakNumber]))
    {
        peaks.erase(peaks.begin() + peakNumber);
        delete gaus;
        return;
    }
}

// Ensure all methods that do not modify the state are marked as const
std::string Histogram::getMainHistName() const
{
    if (!mainHist)
    {
        return "Invalid histogram";
    }
    std::string name = mainHist->GetName();
    std::size_t pos = name.find_last_not_of("0123456789");
    if (pos != std::string::npos && pos + 1 < name.size())
    {
        return name.substr(pos + 1);
    }
    return name;
}
