#include "../include/Histogram.h"

// Variabile globale
float MaxDistance = 10;
float MinDistance = 1.9f;
long int iterations = 0;
long int goodGaus = 0;
long int badGaus = 0;

Histogram::Histogram(int xMin, int xMax, int maxFWHM, float minAmplitude, float maxAmplitude, int numberOfPeaks, TH1D *mainHist, const std::string &TH2histogram_name, std::string sourceName)
    : xMin(xMin), xMax(xMax), maxFWHM(maxFWHM), minAmplitude(minAmplitude), maxAmplitude(maxAmplitude), numberOfPeaks(numberOfPeaks), mainHist(mainHist), tempHist(nullptr), calibratedHist(nullptr), peakCount(0), m(0), b(0), TH2histogram_name(TH2histogram_name), sourceName(sourceName), peakMatchCount(0)
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
    m = histogram.m;
    b = histogram.b;
    peakCount = histogram.peakCount;
    TH2histogram_name = histogram.TH2histogram_name; // Corrected
    sourceName = histogram.sourceName;
    peakMatchCount = histogram.peakMatchCount;
    peaks = histogram.peaks;

    mainHist = (histogram.mainHist) ? new TH1D(*histogram.mainHist) : nullptr;
    tempHist = (histogram.tempHist) ? new TH1D(*histogram.tempHist) : nullptr;
    calibratedHist = (histogram.calibratedHist) ? new TH1D(*histogram.calibratedHist) : nullptr;
}


Histogram &Histogram::operator=(const Histogram &histogram)
{
    if (this != &histogram)
    {
        // Release old resources
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
        m = histogram.m;
        b = histogram.b;
        peakCount = histogram.peakCount;
        TH2histogram_name = histogram.TH2histogram_name; // Corrected
        sourceName = histogram.sourceName;
        peakMatchCount = histogram.peakMatchCount;
        peaks = histogram.peaks;

        // Allocate new resources
        mainHist = (histogram.mainHist) ? new TH1D(*histogram.mainHist) : nullptr;
        tempHist = (histogram.tempHist) ? new TH1D(*histogram.tempHist) : nullptr;
        calibratedHist = (histogram.calibratedHist) ? new TH1D(*histogram.calibratedHist) : nullptr;
    }
    return *this;
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

    // Asigurăm că limitele sunt în interiorul histogramei
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
    peakMatchCount = bestCorrelation;
    refineCalibration(knownEnergies, size);
    polinomDegree = 1;
}
// The function that determinates the degree of the polynomial
// Not done yet

int Histogram::getTheDegreeOfPolynomial() const
{
    if (numberOfPeaks == 0)
    {
        std::cerr << "No peaks available for fitting!" << std::endl;
        return -1;
    }
    std::vector<double> xValues;
    std::vector<double> yValues;

    for (int i = 0; i < numberOfPeaks; i++)
    {
        //std::cout << "Peak number: " << i << std::endl;
        double peakPosition = peaks[i].getPosition();
        double peakAmplitude = peaks[i].getAmplitude();
        xValues.push_back(peakPosition);
        yValues.push_back(peakAmplitude);
        //std::cout << "Peak position: " << peakPosition << ", amplitude: " << peakAmplitude << std::endl;
    }
    TGraph graph(xValues.size(), xValues.data(), yValues.data());

    int bestDegree = 0;
    double bestR2 = -1  * pow(10, 10);
    double lowestChi2 = 1e10;

    for (int degree = 1; degree <= 5; ++degree)
    {
        TF1 fitFunction("fitFunction", ("pol" + std::to_string(degree)).c_str(), xValues.front(), xValues.back());

        graph.Fit(&fitFunction, "Q");

        double chi2 = fitFunction.GetChisquare();
        int ndf = fitFunction.GetNDF();
        double r2 = 1 - (chi2 / (ndf > 0 ? ndf : 1));

        //std::cout << "Degree: " << degree << ", R^2: " << r2 << ", Chi2: " << chi2 << std::endl;

        //std::cout<<"R2: "<<r2<<" BestR2: "<<bestR2<<std::endl;
        //std::cout << (r2 > bestR2) << std::endl;
        //std::cout<<"//////"<<std::endl;
        if (r2 > bestR2)
        {
            bestR2 = r2;
            bestDegree = degree;
            lowestChi2 = chi2;
        }
    }

    //std::cout << "Best polynomial degree: " << bestDegree << " with R^2: " << bestR2 << ", Chi2: " << lowestChi2 << std::endl;

    return bestDegree;
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
            // std::cerr << "Error: Bin index out of range. Bin index: " << bin_original << std::endl;
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

void Histogram::changePeak(int peakNumber, double newPosition)
{
    // Verificăm dacă numărul peak-ului este valid
    if (peakNumber < 0 || peakNumber >= peaks.size())
    {
        std::cerr << "Error: Invalid peak number." << std::endl;
        return;
    }

    // Creăm un fit gaussian cu noua poziție
    TF1 *gaus = createGaussianFit(newPosition);
    if (gaus == nullptr)
    {
        std::cerr << "Error: Failed to create Gaussian fit." << std::endl;
        return;
    }

    // Aplicăm fit-ul pe histogramă
    tempHist->Fit(gaus, "RQ");

    // Eliminăm peak-ul existent
    eliminatePeak(peaks[peakNumber]);

    // Adăugăm noul peak
    peaks[peakNumber] = {gaus, mainHist};
    std::cout << "Peak number: " << peakNumber << std::endl;
    std::cout << "Peak position: " << peaks[peakNumber].getPosition() << std::endl;
    // Verificăm condițiile pentru noul peak
    if (!checkConditions(peaks[peakNumber]))
    {
        // Eliminăm peak-ul dacă condițiile nu sunt îndeplinite
        peaks.erase(peaks.begin() + peakNumber);
        delete gaus; // Curățăm obiectul TF1 pentru a preveni scurgerile de memorie
        return;
    }

    // Dacă peak-ul este valid, nu mai trebuie să ștergem `gaus`
    // deoarece `peaks[peakNumber]` va gestiona viața sa.
}

void Histogram::printHistogramWithPeaksRoot(TFile *outputFile)
{
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return;
    }

    outputFile->cd();

    // Resetăm histogramul pentru a elimina ajustările anterioare
    mainHist->GetListOfFunctions()->Clear();

    // Recalculăm și ajustăm histogramul cu noi funcții gaussiene
    for (int i = 0; i < peaks.size(); ++i)
    {
        double newPosition = peaks[i].getPosition();
        TF1 *gaussianFunction = createGaussianFit(newPosition); // Funcție care creează noua ajustare gaussiană
        if (gaussianFunction)
        {
            // peaks[i].setGaussianFunction(gaussianFunction); // Asumând că există această metodă în clasă

            // Adăugăm funcția de ajustare la histogramă
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

void Histogram::outputPeaksDataJson(std::ofstream &jsonFile)
{
    jsonFile << "{\n";
    jsonFile << "\t\"Source\": \"" << sourceName << "\",\n";
    jsonFile << "\t\"Histogram\": \"" << mainHist->GetName() << "\",\n";
    jsonFile << "\t\"TH2Source_FileName\": \"" << TH2histogram_name << "\",\n";
    jsonFile << "\t\"NumberOfPeaks\": " << peaks.size() << ",\n";
    jsonFile << "\t\"Calibration Degree\": " << getTheDegreeOfPolynomial()/*polinomDegree*/ << ",\n";
    jsonFile << "\t\"Calibration Factor m\": " << m << ",\n";
    jsonFile << "\t\"Calibration Factor b\": " << b << ",\n";
    jsonFile << "\t\"Peaks\": [\n";

    for (size_t i = 0; i < peaks.size(); ++i)
    {
        jsonFile << "\t\t{\n";
        //jsonFile << "\t\t\t\"Number_Peak\": " << i + 1 << ",\n";
        jsonFile << "\t\t\t\"position\": " << peaks[i].getPosition() << ",\n";
        jsonFile << "\t\t\t\"FWHM\": " << peaks[i].getFWHM() << ",\n";
        jsonFile << "\t\t\t\"area\": " << peaks[i].getArea() << "\n";
        // jsonFile << "\t\t\t\"amplitude\": " << peaks[i].getAmplitude() << ",\n";
        // jsonFile << "\t\t\t\"sigma\": " << peaks[i].getSigma() << ",\n";
        // jsonFile << "\t\t\t\"leftLimit\": " << peaks[i].getLeftLimit() << ",\n";
        // jsonFile << "\t\t\t\"rightLimit\": " << peaks[i].getRightLimit() << "\n";

        if (i < peaks.size() - 1)
        {
            jsonFile << "\t\t},\n"; // Virgula dacă nu e ultimul element
        }
        else
        {
            jsonFile << "\t\t}\n"; // Fără virgulă la ultimul element
        }
    }

    jsonFile << "\t]\n";
    jsonFile << "}\n";
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