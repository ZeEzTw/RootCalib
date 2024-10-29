#include "../include/Histogram.h"

// global limits for the peaks
float MaxDistance = 10;
float MinDistance = 1.9f;
long int iterations = 0;
long int goodGaus = 0;
long int badGaus = 0;
Histogram::Histogram()
{
    xMin = 0;
    xMax = 0;
    maxFWHM = 0;
    minAmplitude = 0;
    maxAmplitude = 0;
    numberOfPeaks = 0;
    mainHist = nullptr;
    tempHist = nullptr;
    calibratedHist = nullptr;
    m = 0;
    b = 0;
    polinomDegree = 0;
    peakMatchCount = 0;
    polynomialFitThreshold = 1e-3;
    TH2histogram_name = "An empty histogram";
    sourceName = "Empty source";
    peakCount = 0;
}

Histogram::Histogram(int xMin, int xMax, int maxFWHM, float minAmplitude, float maxAmplitude, std::string serial, int detType, float polynomialFitThreshold, int numberOfPeaks, TH1D *mainHist, const std::string &TH2histogram_name, std::string sourceName)
    : xMin(xMin), xMax(xMax), maxFWHM(maxFWHM), minAmplitude(minAmplitude), maxAmplitude(maxAmplitude), serial(serial), detType(detType), numberOfPeaks(numberOfPeaks), mainHist(mainHist), tempHist(nullptr), calibratedHist(nullptr), peakCount(0), m(0), b(0), TH2histogram_name(TH2histogram_name), sourceName(sourceName), peakMatchCount(0), polynomialFitThreshold(polynomialFitThreshold) // Initialize threshold
{
    if (mainHist != nullptr)
    {
        this->tempHist = (TH1D *)mainHist->Clone();
        this->calibratedHist = (TH1D *)mainHist->Clone();
    }
    else
    {
        std::cerr << "Warning: mainHist is null in Histogram constructor" << std::endl;
    }
}

Histogram::~Histogram()
{
    /*if (tempHist)
    {
        delete tempHist;
    }
    if (calibratedHist)
    {
        delete calibratedHist;
    }*/
}
Histogram::Histogram(const Histogram &histogram)
{
    xMin = histogram.xMin;
    xMax = histogram.xMax;
    maxAmplitude = histogram.maxAmplitude;
    minAmplitude = histogram.minAmplitude;
    maxFWHM = histogram.maxFWHM;
    numberOfPeaks = histogram.numberOfPeaks;
    serial = histogram.serial;
    detType = histogram.detType;
    m = histogram.m;
    b = histogram.b;
    peakCount = histogram.peakCount;
    TH2histogram_name = histogram.TH2histogram_name;
    sourceName = histogram.sourceName;
    peakMatchCount = histogram.peakMatchCount;
    polynomialFitThreshold = histogram.polynomialFitThreshold;

    mainHist = (histogram.mainHist) ? new TH1D(*histogram.mainHist) : nullptr;
    tempHist = (histogram.tempHist) ? new TH1D(*histogram.tempHist) : nullptr;
    calibratedHist = (histogram.calibratedHist) ? new TH1D(*histogram.calibratedHist) : nullptr;
}

Histogram &Histogram::operator=(const Histogram &histogram)
{
    if (this != &histogram)
    {
        // Free old resources
        delete mainHist;
        delete tempHist;
        delete calibratedHist;

        // Copy values
        xMin = histogram.xMin;
        xMax = histogram.xMax;
        maxAmplitude = histogram.maxAmplitude;
        minAmplitude = histogram.minAmplitude;
        maxFWHM = histogram.maxFWHM;
        numberOfPeaks = histogram.numberOfPeaks;
        serial = histogram.serial;
        detType = histogram.detType;
        m = histogram.m;
        b = histogram.b;
        peakCount = histogram.peakCount;
        TH2histogram_name = histogram.TH2histogram_name;
        sourceName = histogram.sourceName;
        peakMatchCount = histogram.peakMatchCount;
        polynomialFitThreshold = histogram.polynomialFitThreshold;

        mainHist = (histogram.mainHist) ? new TH1D(*histogram.mainHist) : nullptr;
        tempHist = (histogram.tempHist) ? new TH1D(*histogram.tempHist) : nullptr;
        calibratedHist = (histogram.calibratedHist) ? new TH1D(*histogram.calibratedHist) : nullptr;
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
    double leftLimit = MinDistance;
    double rightLimit = MinDistance;

    for (int bin = 1; bin <= tempHist->GetNbinsX(); ++bin)
    {
        float binContent = tempHist->GetBinContent(bin);
        if (binContent == 0)
            continue;

        iterations++;
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

    if (peaks.back().getArea() < 0 || (peaks.back().getPosition() == peaks[peaks.size() - 2].getPosition()) || peaks.back().getPosition() < 0)
    {
        /*09.29.24 This check take care of the cases when the fit is so bad, that the eliminationPeak() fail, the real peak is at 12 and the gaus say that is at 150
        the elimination function eliminate the position 150, but the peak is at 12, so the peak is not eliminated, will repeat
        And the area is negative when again the fit is so bad that the peak is not a peak, is a valley, so the area is negative
        */
        peaks.pop_back();
        delete gaus;
        return -1;
    }
    if (!checkConditions(peaks.back()))
    {
        eliminatePeak(peaks.back());
        peaks.pop_back();
        delete gaus;
        return 0;
    }

    peakCount++;
    eliminatePeak(peaks.back());

    delete gaus;
    return 0;
}

double Histogram::refineCalibration(const double knownEnergies[], int size) const
{
    double totalError = 0.0;

    for (const auto &peak : peaks)
    {
        double predictedEnergy = m * peak.getPosition() + b;
        double minError = std::numeric_limits<double>::max();

        for (int i = 0; i < size; ++i)
        {
            double error = std::abs(predictedEnergy - knownEnergies[i]);
            minError = std::min(minError, error);
        }

        totalError += minError;
    }

    return totalError / peaks.size();
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
    /*for(int i = 0; i < peaks.size(); i++)
    {
        std::cout << "Peak " << peaks[i].getPosition() << " associated with " << peaks[i].getAssociatedPosition() << std::endl;
    }
    */
    peakMatchCount = bestCorrelation;
    calibratePeaksByDegree(); // mai bun cu acesta
    // getTheDegreeOfPolynomial();
    //   Continuă cu calibrarea
    //   getTheDegreeOfPolynomial();
    //   calibratePeaksByDegree();
}

double Histogram::refineCalibrationM()
{
    double mediumM = 0.0;
    int peaksRemained = 0;
    for (const auto &peak : peaks)
    {
        if (peak.getAssociatedPosition() < 1)
        {
            continue;
        }
        mediumM += peak.getAssociatedPosition() / peak.getPosition();
        peaksRemained++;
    }
    return mediumM / peaksRemained;
}

double Histogram::refineCalibrationB()
{
    double totalError = 0.0;
    int peaksRemained = 0;
    for (const auto &peak : peaks)
    {
        double predictedEnergy = m * peak.getPosition() + b;
        if (peak.getAssociatedPosition() == 0.0)
        {
            continue;
        }
        double error = std::abs(predictedEnergy - peak.getAssociatedPosition());
        totalError += error;
        peaksRemained++;
    }

    return totalError / peaksRemained;
}

// getting polynomial degree + values
#include <vector>
#include <cmath>
#include <iostream>

void Histogram::calibratePeaksByDegree()
{
    degree = 1;

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
        std::cerr << "No valid peaks for calibration!" << std::endl;
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

        std::vector<std::vector<double>> XtX = multiplyTransposeMatrix(X);
        std::vector<double> XtY = multiplyTransposeVector(X, Y);
        std::vector<double> coeffs = solveSystem(XtX, XtY);

        if (std::abs(coeffs[currentDegree]) >= polynomialFitThreshold)
        {
            degree = currentDegree;
            coefficients = coeffs;
        }
    }
}

// Helper function to multiply a matrix by its transpose
std::vector<std::vector<double>> Histogram::multiplyTransposeMatrix(const std::vector<std::vector<double>> &X)
{
    int rows = X.size();
    int cols = X[0].size();
    std::vector<std::vector<double>> XtX(cols, std::vector<double>(cols, 0.0));

    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < rows; ++k)
            {
                sum += X[k][i] * X[k][j];
            }
            XtX[i][j] = XtX[j][i] = sum;
        }
    }

    return XtX;
}

// Helper function to multiply the transpose of a matrix by a vector
std::vector<double> Histogram::multiplyTransposeVector(const std::vector<std::vector<double>> &X, const std::vector<double> &Y)
{
    int rows = X.size();
    int cols = X[0].size();
    std::vector<double> XtY(cols, 0.0);

    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j < rows; ++j)
        {
            XtY[i] += X[j][i] * Y[j];
        }
    }

    return XtY;
}

// Helper function to solve a linear system using Gaussian elimination
std::vector<double> Histogram::solveSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{
    int n = A.size();
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(n + 1));

    // Augment the matrix A with the vector b
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = b[i];
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i)
    {
        // Partial pivoting (find the row with the largest element in the current column)
        double maxElement = std::abs(augmentedMatrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
        {
            if (std::abs(augmentedMatrix[k][i]) > maxElement)
            {
                maxElement = std::abs(augmentedMatrix[k][i]);
                maxRow = k;
            }
        }

        // Swap rows if necessary
        if (maxRow != i)
        {
            std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);
        }

        // Make the diagonal element 1 and eliminate below
        for (int k = i + 1; k < n; ++k)
        {
            double c = -augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j <= n; ++j)
            {
                if (i == j)
                {
                    augmentedMatrix[k][j] = 0;
                }
                else
                {
                    augmentedMatrix[k][j] += c * augmentedMatrix[i][j];
                }
            }
        }
    }

    // Back substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; --k)
        {
            augmentedMatrix[k][n] -= augmentedMatrix[k][i] * x[i];
        }
    }

    return x;
}

// V2 calibration section

// Funție pentru validarea și extragerea datelor pentru fit
bool Histogram::extractDataForFit(std::vector<double> &xValues, std::vector<double> &yValues)
{
    if (numberOfPeaks == 0)
    {
        std::cerr << "No peaks available for fitting!" << std::endl;
        degree = 0;
        return false;
    }
    if (numberOfPeaks < 2)
    {
        std::cerr << "Insufficient points for fitting a polynomial." << std::endl;
        degree = 0;
        return false;
    }

    for (const auto &peak : peaks)
    {
        if (peak.getAssociatedPosition() != 0)
        {
            xValues.push_back(peak.getPosition());
            yValues.push_back(peak.getAssociatedPosition());
        }
    }

    if (xValues.size() < 2)
    {
        std::cerr << "Not enough valid peaks for polynomial fitting." << std::endl;
        return false;
    }

    return true;
}

bool Histogram::areCoefficientsValid(TF1 *fitFunction, int degree, double threshold)
{
    for (int i = degree; i >= 0; --i)
    {
        double coefficient = fitFunction->GetParameter(i);
        if (std::abs(coefficient) < threshold)
        {
            std::cout << "Coefficient a" << i << " = " << fitFunction->GetParameter(i) << " is below threshold." << std::endl;
            return false;
        }
    }
    return true;
}

void Histogram::getTheDegreeOfPolynomial()
{
    std::vector<double> xValues, yValues;
    if (!extractDataForFit(xValues, yValues))
        return;

    TGraph graph(xValues.size(), xValues.data(), yValues.data());
    graph.SetTitle("Polynomial Fit;X-axis;Y-axis");
    graph.SetMarkerStyle(20);
    graph.SetMarkerColor(kBlue);

    TFile outFile("polynomial_fits.root", "UPDATE");

    int degree = 1; // Start with degree 1
    coefficients.clear();
    TF1 *fitFunction = nullptr;

    bool foundValidFit = false; // Flag to check if a valid fit is found

    while (degree <= std::min(3, numberOfPeaks - 1))
    {
        fitFunction = new TF1(("fitFunc_deg" + std::to_string(degree)).c_str(),
                              ("pol" + std::to_string(degree)).c_str(),
                              xValues.front(), xValues.back());

        if (graph.Fit(fitFunction, "RQ") != 0)
        {
            std::cerr << "Fit failed for degree: " << degree << std::endl;
            delete fitFunction;
            break;
        }

        fitFunction->Write(("PolyFit_deg" + std::to_string(degree)).c_str());

        // Check the validity of the coefficients
        if (areCoefficientsValid(fitFunction, degree, polynomialFitThreshold))
        {
            coefficients.clear(); // If valid, update coefficients
            for (int i = 0; i <= degree; ++i)
            {
                coefficients.push_back(fitFunction->GetParameter(i));
            }
            foundValidFit = true; // Mark that a valid fit was found
            ++degree;             // Increase the degree for the next check
        }
        else
        {
            std::cout << "Stopping at degree " << degree - 1 << " due to small coefficients." << std::endl;
            delete fitFunction;
            break;
        }

        delete fitFunction;
    }

    // If no valid fit is found, use a default polynomial of degree 1
    if (!foundValidFit && coefficients.empty())
    {
        std::cout << "Using polynomial of degree 1 as the default fit." << std::endl;
        fitFunction = new TF1("fitFunc_deg1", "pol1", xValues.front(), xValues.back());
        if (graph.Fit(fitFunction, "RQ") == 0)
        {
            for (int i = 0; i <= 1; ++i)
            {
                coefficients.push_back(fitFunction->GetParameter(i));
            }
        }
        delete fitFunction;
    }

    graph.Write(("GraphWithDataPoints_" + getMainHistName()).c_str());
    outFile.Close();

    std::cout << "Best polynomial degree: " << (coefficients.size() - 1) << std::endl;
    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        std::cout << "Coefficient a" << i << " = " << coefficients[i] << std::endl;
    }
}

// Apply calibration section
void Histogram::initializeCalibratedHist()
{
    if (mainHist == nullptr)
    {
        std::cerr << "Error: mainHist is not initialized." << std::endl;
        return;
    }

    if (calibratedHist == nullptr)
    {
        int numBins = mainHist->GetNbinsX();
        calibratedHist = new TH1D("calibratedHist", "Calibrated Histogram", numBins, mainHist->GetXaxis()->GetXmin(), mainHist->GetXaxis()->GetXmax());
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
        jsonFile << "\t\t\t\t\"area\": [" << peaks[i].getArea() << ", " << peaks[1].getAreaError() << "]\n";

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
}

void Histogram::printHistogramWithPeaksRoot(TFile *outputFile)
{
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: Could not open file for writing" << std::endl;
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
        std::cerr << "Error: Could not open file for writing" << std::endl;
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

    if (left > MaxDistance || left < MinDistance)
    {
        badGaus++;
        goodGaus--;
        leftLimitPosition = peak.getPosition() - MinDistance;
    }
    if (right > MaxDistance || right < MinDistance)
    {
        rightLimitPosition = peak.getPosition() + MinDistance;
    }

    goodGaus++;
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

unsigned int Histogram::getpeakMatchCount() const
{
    return peakMatchCount;
}

TH1D *Histogram::getCalibratedHist() const
{
    return calibratedHist;
}

TH1D *Histogram::getMainHist() const
{
    return mainHist;
}

std::string Histogram::getMainHistName() const
{
    std::string name = mainHist->GetName();

    std::size_t pos = name.find_last_not_of("0123456789");

    if (pos != std::string::npos && pos + 1 < name.size())
    {
        return name.substr(pos + 1);
    }

    return mainHist->GetName();
}

// extra function section
void Histogram::changePeak(int peakNumber, double newPosition)
{
    {
        std::cerr << "Error: Invalid peak number." << std::endl;
        return;
    }

    TF1 *gaus = createGaussianFit(newPosition);
    if (gaus == nullptr)
    {
        std::cerr << "Error: Failed to create Gaussian fit." << std::endl;
        return;
    }

    tempHist->Fit(gaus, "RQ");

    eliminatePeak(peaks[peakNumber]);

    peaks[peakNumber] = {gaus, mainHist};
    std::cout << "Peak number: " << peakNumber << std::endl;
    std::cout << "Peak position: " << peaks[peakNumber].getPosition() << std::endl;
    if (!checkConditions(peaks[peakNumber]))
    {
        peaks.erase(peaks.begin() + peakNumber);
        delete gaus;
        return;
    }
}

/*void Histogram::initializeCalibratedHist()
{
    if (mainHist == nullptr)
    {
        std::cerr << "Error: mainHist is not initialized." << std::endl;
        return;
    }
    if (calibratedHist == nullptr)
    {
        int numBins = mainHist->GetNbinsX();
        calibratedHist = new TH1D("calibratedHist", "Calibrated Histogram", numBins, mainHist->GetXaxis()->GetXmin(), mainHist->GetXaxis()->GetXmax());
    }
}
double Histogram::getInterpolatedContent(int bin_original, double binCenter_original) const
{
    double content_original = mainHist->GetBinContent(bin_original);
    double content_residual = (bin_original < mainHist->GetNbinsX()) ? mainHist->GetBinContent(bin_original + 1) : content_original;
    double binLowEdge_original = mainHist->GetXaxis()->GetBinLowEdge(bin_original);
    double binWidth_original = mainHist->GetXaxis()->GetBinWidth(bin_original);
    double fraction = (binCenter_original - binLowEdge_original) / binWidth_original;
    return content_original + fraction * (content_residual - content_original);
}

double Histogram::inversePolynomial(double y) const
{
    // Pentru gradul 1, putem calcula direct
    if (degree == 1)
    {
        return (y - coefficients[0]) / coefficients[1];
    }

    // Pentru grade superioare, folosim căutare binară
    double xMin = mainHist->GetXaxis()->GetXmin();
    double xMax = mainHist->GetXaxis()->GetXmax();
    double epsilon = 1e-6; // Precizia dorită

    while (xMax - xMin > epsilon)
    {
        double xMid = (xMin + xMax) / 2;
        double yMid = evaluatePolynomial(xMid);

        if (std::abs(yMid - y) < epsilon)
        {
            return xMid;
        }

        if (yMid < y)
        {
            xMin = xMid;
        }
        else
        {
            xMax = xMid;
        }
    }

    return (xMin + xMax) / 2;
}

// Funcția evaluatePolynomial rămâne neschimbată
double Histogram::evaluatePolynomial(double x) const
{
    double result = 0.0;
    for (int i = 0; i <= degree; ++i)
    {
        result += coefficients[i] * std::pow(x, i);
    }
    return result;
}

void Histogram::applyXCalibration()
{
    initializeCalibratedHist();
    int numBinsCalibrated = calibratedHist->GetNbinsX();
    int numBinsOriginal = mainHist->GetNbinsX();

    for (int bin_calibrated = 1; bin_calibrated <= numBinsCalibrated; ++bin_calibrated)
    {
        double binCenter_calibrated = calibratedHist->GetXaxis()->GetBinCenter(bin_calibrated);
        double binCenter_original = inversePolynomial(binCenter_calibrated);

        int bin_original = mainHist->GetXaxis()->FindBin(binCenter_original);
        if (bin_original < 1 || bin_original >= numBinsOriginal)
        {
            continue;
        }

        double interpolated_content = getInterpolatedContent(bin_original, binCenter_original);
        calibratedHist->SetBinContent(bin_calibrated, interpolated_content);
    }
}

void Histogram::applyXCalibration()
{
    if (mainHist == nullptr)
    {
        std::cerr << "Error: mainHist is not initialized." << std::endl;
        return;
    }
    if (calibratedHist == nullptr)
    {
        int numBins = mainHist->GetNbinsX();
        calibratedHist = new TH1D("calibratedHist", "Calibrated Histogram", numBins, mainHist->GetXaxis()->GetXmin(), mainHist->GetXaxis()->GetXmax());
    }

    int numBinsCalibrated = calibratedHist->GetNbinsX();
    int numBinsOriginal = mainHist->GetNbinsX();

    for (int bin_calibrated = 1; bin_calibrated <= numBinsCalibrated; ++bin_calibrated)
    {
        double binCenter_calibrated = calibratedHist->GetXaxis()->GetBinCenter(bin_calibrated);
        double binCenter_original = (binCenter_calibrated - b) / m;
        int bin_original = mainHist->GetXaxis()->FindBin(binCenter_original);
        if (bin_original < 1 || bin_original >= numBinsOriginal)
        {
            std::cerr << "Error: Bin index out of range." << std::endl;
            continue;
        }

        double content_original = mainHist->GetBinContent(bin_original);
        double content_residual = (bin_original < numBinsOriginal) ? mainHist->GetBinContent(bin_original + 1) : content_original;
        double binLowEdge_original = mainHist->GetXaxis()->GetBinLowEdge(bin_original);
        double binWidth_original = mainHist->GetXaxis()->GetBinWidth(bin_original);
        double fraction = (binCenter_original - binLowEdge_original) / binWidth_original;
        double interpolated_content = content_original + fraction * (content_residual - content_original);
        calibratedHist->SetBinContent(bin_calibrated, interpolated_content);
    }
}*/

/*void Histogram::calibratePeaks(const double knownEnergies[], int size)
{
    double bestM = 0.0;
    double bestB = 0.0;
    int bestCorrelation = 0;
    double valueAssociatedWith = 0.0;

    for (double m = 1.0; m <= 5.0; m += 0.0001)
    {
        std::vector<double> associatedValues(peaks.size(), 0.0); // Inițializează vectorul cu dimensiunea potrivită
        int correlations = 0;
        int peakCount = 0;

        for (const auto &peak : peaks)
        {
            double predictedEnergy = m * peak.getPosition() + b;
            if (checkPredictedEnergies(predictedEnergy, knownEnergies, size, 5, valueAssociatedWith))
            {
                ++correlations;
                associatedValues[peakCount] = valueAssociatedWith;
            }
            else
            {
                // associatedValues[peakCount] = 0.0;
            }
            peakCount++;
        }

        if (correlations > bestCorrelation)
        {
            bestCorrelation = correlations;
            int i = 0;
            for (auto &peak : peaks)
            {
                peak.setAssociatedPosition(associatedValues[i]);
                ++i;
            }

            this->m = m;
        }
    }

    peakMatchCount = bestCorrelation;
    m = bestM;
    // b = bestB;
    // m = refineCalibrationM();
    b = refineCalibrationB();
    getTheDegreeOfPolynomial();
    calibratePeaksByDegree();
    polinomDegree = 1;
}*/

/*std::cout << "Coeficienții polinomului sunt: ";
for (int i = 0; i <= degree; ++i)
{
    std::cout << "a" << i << " = " << coeffs(i) << ", ";
}
std::cout << std::endl;

for (auto &peak : peaks)
{
    if (peak.getAssociatedPosition() != 0)
    {
        double x = peak.getPosition();
        if (degree == 2)
        {
            double calibratedEnergy = coeffs(2) * x * x + coeffs(1) * x + coeffs(0);
            std::cout << "Pozitie: " << x << ", Energie calibrata: " << calibratedEnergy << std::endl;
        }
        else if (degree == 3)
        {
            double calibratedEnergy = coeffs(3) * x * x * x + coeffs(2) * x * x + coeffs(1) * x + coeffs(0);
            std::cout << "Pozitie: " << x << ", Energie calibrata: " << calibratedEnergy << std::endl;
        }
        else if (degree == 4)
        {
            double calibratedEnergy = coeffs(4) * x * x * x * x + coeffs(3) * x * x * x + coeffs(2) * x * x + coeffs(1) * x + coeffs(0);
            std::cout << "Pozitie: " << x << ", Energie calibrata: " << calibratedEnergy << std::endl;
        }
        else
        {
            double calibratedEnergy = coeffs(1) * x + coeffs(0);
        }
    }
}

}

//applyXCalibration old version

double newtonRaphson(const Eigen::VectorXd &coeffs, double initialGuess, int maxIter = 100, double tol = 1e-6)
{
    double x = initialGuess;

    for (int iter = 0; iter < maxIter; ++iter)
    {
        // Evaluate the polynomial and its derivative
        double fx = 0.0;
        double fpx = 0.0;

        int degree = coeffs.size() - 1;

        for (int i = 0; i <= degree; ++i)
        {
            fx += coeffs[i] * std::pow(x, degree - i);
            if (i < degree)
            {
                fpx += (degree - i) * coeffs[i] * std::pow(x, degree - i - 1);
            }
        }

        // Check if the derivative is too small (to avoid division by zero)
        if (std::abs(fpx) < tol)
        {
            std::cerr << "Derivative too small during Newton-Raphson iteration." << std::endl;
            break;
        }

        // Newton-Raphson step
        double x_new = x - fx / fpx;

        // Check for convergence
        if (std::abs(x_new - x) < tol)
        {
            return x_new;
        }

        x = x_new;
    }

    std::cerr << "Newton-Raphson method did not converge." << std::endl;
    return x; // Return the last estimate
}
// Trebuie sa ma ocup aici nu este bine sigur ce fac
double Histogram::getInterpolatedContent(int bin_original, double binCenter_original) const
{
    double content_original = mainHist->GetBinContent(bin_original);
    double content_residual = (bin_original < mainHist->GetNbinsX()) ? mainHist->GetBinContent(bin_original + 1) : content_original;

    double binLowEdge_original = mainHist->GetXaxis()->GetBinLowEdge(bin_original);
    double binWidth_original = mainHist->GetXaxis()->GetBinWidth(bin_original);
    double fraction = (binCenter_original - binLowEdge_original) / binWidth_original;

    return content_original + fraction * (content_residual - content_original);
}
*/