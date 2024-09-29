#include "../include/Histogram.h"

// Variabile globale
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
    TH2histogram_name = "An empty histogram";
    sourceName = "Empty source";
    peakCount = 0;
}

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
    int result = 0;
    int count = 0;
    std::cout<<"Number of peaks: "<<numberOfPeaks<<std::endl;
    while (result == 0 && count < numberOfPeaks) // Căutăm vârfuri până găsim un vârf invalid sau ajungem la 20 de încercări.
    {
        result = detectAndFitPeaks();
        if (result == -1)
        {
            std::cout << "BREAKKKKKKKKK" << std::endl;
            break;
        }
        count++;
    }
    // for (int i = 0; i < numberOfPeaks; i++)
    //{
    // detectAndFitPeaks();
    //}
}
void Histogram::eliminatePeak(const Peak &peak)
{
    double mean = peak.getMean();
    double sigma = peak.getSigma();
    int leftLimit = static_cast<int>(mean - 3 * sigma);
    int rightLimit = static_cast<int>(mean + 3 * sigma);

    // Asigurăm că limitele sunt în interiorul histogramei
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
        /*09.29.24 This check take care of the cases when the fit is so bad, that the eliminationPeak() fail, tha real peak is at 12 and the gaus say that is at 150
        the elimination function eliminate the position 150, but the peak is at 12, so the peak is not eliminated, will repeat
        And the are is negative when again the fit is so bad that the peak is not a peak, is a valley, so the area is negative
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

    for (double m = 1.0; m <= 5.0; m += 0.0001)
    {
        for (double b = 0.0; b <= 0.0; b += 0.1)
        {
            std::vector<double> associatedValues(peaks.size(), 0.0); // Inițializează vectorul cu dimensiunea potrivită
            int correlations = 0;
            int peakCount = 0;

            for (const auto &peak : peaks)
            {
                double predictedEnergy = m * peak.getPosition() + b;
                if (checkPredictedEnergies(predictedEnergy, knownEnergies, size, 2, valueAssociatedWith))
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
                this->b = b;
            }
        }
    }

    peakMatchCount = bestCorrelation;
    // m = bestM;
    // b = bestB;
    // m = refineCalibrationM();
    b = refineCalibrationB();
    polinomDegree = 1;
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
        std::cout << peak.getAssociatedPosition() << " " << peak.getPosition() << " " << (peak.getAssociatedPosition() / peak.getPosition()) << std::endl;
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
        // std::cout << "Peak number: " << i << std::endl;
        double peakPosition = peaks[i].getPosition();
        double peakAmplitude = peaks[i].getAmplitude();
        xValues.push_back(peakPosition);
        yValues.push_back(peakAmplitude);
        // std::cout << "Peak position: " << peakPosition << ", amplitude: " << peakAmplitude << std::endl;
    }
    TGraph graph(xValues.size(), xValues.data(), yValues.data());

    int bestDegree = 0;
    double bestR2 = -1 * pow(10, 10);
    double lowestChi2 = 1e10;

    for (int degree = 1; degree <= 5; ++degree)
    {
        TF1 fitFunction("fitFunction", ("pol" + std::to_string(degree)).c_str(), xValues.front(), xValues.back());

        graph.Fit(&fitFunction, "Q");

        double chi2 = fitFunction.GetChisquare();
        int ndf = fitFunction.GetNDF();
        double r2 = 1 - (chi2 / (ndf > 0 ? ndf : 1));

        // std::cout << "Degree: " << degree << ", R^2: " << r2 << ", Chi2: " << chi2 << std::endl;

        // std::cout<<"R2: "<<r2<<" BestR2: "<<bestR2<<std::endl;
        // std::cout << (r2 > bestR2) << std::endl;
        // std::cout<<"//////"<<std::endl;
        if (r2 > bestR2)
        {
            bestR2 = r2;
            bestDegree = degree;
            lowestChi2 = chi2;
        }
    }

    // std::cout << "Best polynomial degree: " << bestDegree << " with R^2: " << bestR2 << ", Chi2: " << lowestChi2 << std::endl;

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
            if (peaks[i].getAssociatedPosition() == 0)
            {
                gaussianFunction->SetLineColor(kGreen);
            }
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
    jsonFile << "\t\"Histogram\": \"" << getMainHistName() << "\",\n";
    jsonFile << "\t\"TH2: \": \"" << TH2histogram_name << "\",\n";
    jsonFile << "\t\"NumberOfPeaks\": " << peaks.size() << ",\n";
    jsonFile << "\t\"Calibration Degree\": " << getTheDegreeOfPolynomial() /*polinomDegree*/ << ",\n";
    std::vector<float> calibrationFactors = {m, b};

    jsonFile << "\t\"Calibration Data\": [\n"; // Corectat pentru a începe vectorul de date de calibrare

    for (size_t i = 0; i < calibrationFactors.size(); ++i)
    {
        jsonFile << "\t\t" << calibrationFactors[i];
        if (i < calibrationFactors.size() - 1) // Nu adăuga o virgulă după ultimul element
        {
            jsonFile << ",";
        }
        jsonFile << "\n";
    }

    jsonFile << "\t],\n"; // Încheie vectorul de date de calibrare
    jsonFile << "\t\"Peaks\": [\n";

    for (size_t i = 0; i < peaks.size(); ++i)
    {
        jsonFile << "\t\t{\n";
        // jsonFile << "\t\t\t\"Number_Peak\": " << i + 1 << ",\n";
        jsonFile << "\t\t\t\"position\": " << peaks[i].getPosition() << ",\n";
        jsonFile << "\t\t\t\"AssociatedPosition\": " << peaks[i].getAssociatedPosition() << ",\n";
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
