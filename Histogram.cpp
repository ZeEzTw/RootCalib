#include "Histogram.h"

// Variabile globale
float MaxDistance = 10;
float MinDistance = 1.9f;
long int iterations = 0;
long int goodGaus = 0;
long int badGaus = 0;

Histogram::Histogram(int xMin, int xMax, int maxFWHM, int numberOfPeaks, TH1D *mainHist)
    : xMin(xMin), xMax(xMax), maxFWHM(maxFWHM), numberOfPeaks(numberOfPeaks), mainHist(mainHist), tempHist(nullptr), calibratedHist(nullptr), peakCount(0), m(0), b(0)
{

    if (mainHist != nullptr)
    {
        this->tempHist = (TH1D *)mainHist->Clone(); // Clonăm histogramă pentru procesare ulterioară
    }
    else
    {
        std::cerr << "Warning: mainHist is null in Histogram constructor" << std::endl;
    }
}

Histogram::~Histogram()
{
    if (tempHist != nullptr)
    {
        delete tempHist;
    }
    if (calibratedHist != nullptr)
    {
        delete calibratedHist;
    }
}
void Histogram::findPeaks()
{
    for (int i = 0; i < numberOfPeaks; i++)
    {
        detectAndFitPeaks();
    }
}
void Histogram::eliminatePeak(const Peak &peak)
{
    double mean = peak.getMean();
    double sigma = peak.getSigma();
    int leftLimit = static_cast<int>(mean - 3 * sigma);
    int rightLimit = static_cast<int>(mean + 3 * sigma);

    // Asigurăm că limitele sunt în interiorul histogramului
    if (leftLimit < 1)
        leftLimit = 1;
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

    return condition1 && condition2;
}
TF1 *Histogram::createGaussianFit(int maxBin)
{
    float maxPeakX = mainHist->GetXaxis()->GetBinCenter(maxBin);
    // is good to use [0]*exp(-0.5*((x-[1])/[2])**2) + [3] + ([4]*x) to fit the gaussian
    TF1 *gaus = new TF1(Form("gausFit_%d", peakCount), "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + ([4]*x) + ([5]*x*x)", maxPeakX - 10, maxPeakX + 10);
    gaus->SetParameters(tempHist->GetBinContent(maxBin), maxPeakX, 1.0, 0.0, 0.0, 0.0); // Initial parameter adjustments
    gaus->SetParLimits(2, 0.1, 10.0);                                                   // Limits for standard deviation
    return gaus;
}

int Histogram::findMaxBin()
{
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
        double peakWithoutBackgroundTemp = binContent - (mainHist->GetBinContent(mainHist->FindBin(bin - leftLimit)) + mainHist->GetBinContent(mainHist->FindBin(bin + rightLimit))) / 2;

        if (peakWithoutBackgroundTemp > peakWithoutBackground)
        {
            peakWithoutBackground = peakWithoutBackgroundTemp;
            maxPeakY = binContent;
            maxBin = bin;
        }
    }

    return maxBin;
}

void Histogram::detectAndFitPeaks()
{
    int maxBin = findMaxBin();
    if (maxBin == 0)
    {
        return;
    }
    TF1 *gaus = createGaussianFit(maxBin);
    tempHist->Fit(gaus, "RQ");
    peaks.emplace_back(gaus, mainHist);
    if (!checkConditions(peaks.back()))
    {
        peaks.pop_back();
    }
    else
    {
        eliminatePeak(peaks.back());
        peakCount++;
    }
    delete gaus;
}

bool Histogram::checkPredictedEnergies(double predictedEnergy, const double knownEnergies[], int size, float errorAdmitted) const
{
    double minError = std::numeric_limits<double>::max();

    for (int i = 0; i < size; ++i)
    {
        double error = std::abs(predictedEnergy - knownEnergies[i]);
        minError = std::min(minError, error);
    }

    return minError < errorAdmitted;
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

void Histogram::calibratePeaks(const double knownEnergies[], int size)
{
    double bestM = 0.0;
    double bestB = 0.0;
    int bestCorrelation = 0;

    for (double m = 0.01; m <= 5.0; m += 0.01)
    {
        for (double b = 0.0; b <= 15.0; b += 0.1)
        {
            int correlations = 0;
            for (const auto &peak : peaks)
            {
                double predictedEnergy = m * peak.getPosition() + b;
                if (checkPredictedEnergies(predictedEnergy, knownEnergies, size, 3))
                {
                    ++correlations;
                }
            }

            if (correlations > bestCorrelation)
            {
                bestCorrelation = correlations;
                this->m = m;
                this->b = b;
            }
        }
    }
    refineCalibration(knownEnergies, size);
}
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
double Histogram::getInterpolatedContent(int bin_original, double binCenter_original) const
{
    double content_original = mainHist->GetBinContent(bin_original);
    double content_residual = (bin_original < mainHist->GetNbinsX()) ? mainHist->GetBinContent(bin_original + 1) : content_original;

    double binLowEdge_original = mainHist->GetXaxis()->GetBinLowEdge(bin_original);
    double binWidth_original = mainHist->GetXaxis()->GetBinWidth(bin_original);
    double fraction = (binCenter_original - binLowEdge_original) / binWidth_original;

    return content_original + fraction * (content_residual - content_original);
}
void Histogram::applyXCalibration()
{
    initializeCalibratedHist();

    int numBinsCalibrated = calibratedHist->GetNbinsX();
    int numBinsOriginal = mainHist->GetNbinsX();

    for (int bin_calibrated = 1; bin_calibrated <= numBinsCalibrated; ++bin_calibrated)
    {
        double binCenter_calibrated = calibratedHist->GetXaxis()->GetBinCenter(bin_calibrated);
        double binCenter_original = (binCenter_calibrated - b) / m;
        int bin_original = mainHist->GetXaxis()->FindBin(binCenter_original);

        if (bin_original < 1 || bin_original >= numBinsOriginal)
        {
            std::cerr << "Error: Bin index out of range. Bin index: " << bin_original << std::endl;
            continue;
        }

        double interpolated_content = getInterpolatedContent(bin_original, binCenter_original);
        calibratedHist->SetBinContent(bin_calibrated, interpolated_content);
    }
}

/*void Histogram::applyXCalibration()
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

void Histogram::printHistogramWithPeaksRoot(TFile *outputFile) const
{
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return;
    }
    outputFile->cd();
    for (int i = 0; i < peaks.size(); ++i)
    {
        TF1 *gaussianFunction = peaks[i].getGaussianFunction();
        if (gaussianFunction)
        {
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

void Histogram::outputPeaksDataJson(std::ofstream &file)
{
    if (!file.is_open())
    {
        std::cerr << "Error: File not open for writing" << std::endl;
        return;
    }

    file << "{\n";
    file << "\t\"peaks\": [\n";

    for (size_t i = 0; i < peaks.size(); ++i)
    {
        peaks[i].outputDataJson(file);
        if (i < peaks.size() - 1)
        {
            file << ",\n";
        }
    }

    file << "\n\t]\n";
    file << "}\n";
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