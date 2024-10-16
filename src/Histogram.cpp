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
    polynomialFitThreshold = 1e-5;
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

    // Iterație prin valori posibile pentru m și b
    for (double m = 0.01; m <= 5.0; m += 0.0001)
    {
        std::vector<double> associatedValues(peaks.size(), 0.0);
        int correlations = 0;
        int peakCount = 0;

        // Iterează prin fiecare peak
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
    degree = 1;
    peakMatchCount = bestCorrelation;
    calibratePeaksByDegree(); // mai bun cu acesta
    // Continuă cu calibrarea
    // getTheDegreeOfPolynomial();
    // calibratePeaksByDegree();
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

void Histogram::getTheDegreeOfPolynomial()
{
    if (numberOfPeaks == 0)
    {
        std::cerr << "No peaks available for fitting!" << std::endl;
        degree = 0;
        return;
    }
    else if (numberOfPeaks < 2)
    {
        degree = 0;
        std::cerr << "Insufficient points for fitting a polynomial." << std::endl;
        return;
    }

    // Inițializează datele pentru ajustare
    std::vector<double> xValues;
    std::vector<double> yValues;

    for (int i = 0; i < numberOfPeaks; i++)
    {
        if (peaks[i].getAssociatedPosition() == 0)
        {
            continue;
        }
        double peakPosition = peaks[i].getPosition();
        double peakAmplitude = peaks[i].getAssociatedPosition();
        xValues.push_back(peakPosition);
        yValues.push_back(peakAmplitude);
    }

    if (xValues.size() < 2)
    {
        std::cerr << "Not enough valid peaks for polynomial fitting." << std::endl;
        return;
    }

    TGraph graph(xValues.size(), xValues.data(), yValues.data());
    graph.SetTitle("Polynomial Fit;X-axis;Y-axis");
    graph.SetMarkerStyle(20);
    graph.SetMarkerColor(kBlue);

    TFile outFile("polynomial_fits.root", "UPDATE");

    int bestDegree = 0;
    double bestR2 = -1e10;
    double lowestChi2 = 1e10;
    TF1 *bestFitFunction = nullptr;

    const double epsilon = 1e-20;

    for (int degree = 1; degree <= std::min(2, numberOfPeaks - 1); ++degree)
    {
        TF1 *fitFunction = new TF1(("fitFunction_" + getMainHistName() + "_deg" + std::to_string(degree)).c_str(),
                                   ("pol" + std::to_string(degree)).c_str(), xValues.front(), xValues.back());

        int fitStatus = graph.Fit(fitFunction, "RQ");
        if (fitStatus != 0)
        {
            std::cerr << "Fit failed for degree: " << degree << std::endl;
            delete fitFunction;
            continue;
        }

        fitFunction->Write(("PolynomialFitDegree_" + std::to_string(degree) + "_" + getMainHistName()).c_str());

        double chi2 = fitFunction->GetChisquare();
        int ndf = fitFunction->GetNDF();
        double r2 = 1 - (chi2 / (ndf > 0 ? ndf : 1));

        bool coefficientsTooSmall = false;
        for (int i = degree; i >= 0; --i)
        {
            double coeff = fitFunction->GetParameter(i);
            if (std::abs(coeff) < epsilon)
            {
                coefficientsTooSmall = true;
                std::cout << "Coefficient a" << i << " is below threshold: " << coeff << std::endl;
                break;
            }
        }

        if (coefficientsTooSmall)
        {
            std::cout << "Polynomial degree " << degree << " rejected due to small coefficients." << std::endl;
            delete fitFunction;
            break;
        }

        if (r2 > bestR2)
        {
            bestR2 = r2;
            bestDegree = degree;
            lowestChi2 = chi2;

            if (bestFitFunction != nullptr)
                delete bestFitFunction;
            bestFitFunction = fitFunction;
        }
        else
        {
            delete fitFunction;
        }
    }

    graph.Write(("GraphWithDataPoints_" + getMainHistName()).c_str());
    outFile.Close();

    if (bestFitFunction != nullptr)
    {
        std::cout << "Best polynomial degree: " << bestDegree << std::endl;
        degree = bestDegree;
        std::cout << "Chi-squared: " << lowestChi2 << std::endl;

        std::cout << "Calibration parameters for the best fit:" << std::endl;
        for (int i = 0; i <= bestDegree; ++i)
        {
            std::cout << "Coefficient a" << i << " = " << bestFitFunction->GetParameter(i) << std::endl;
            coefficients.push_back(bestFitFunction->GetParameter(i));
        }

        delete bestFitFunction;
    }
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
        double binCenter_original;

        if (degree == 1)
        {
            binCenter_original = (binCenter_calibrated - coefficients[0]) / coefficients[1];
        }
        else if (degree == 2)
        {
            double a = coefficients[2];
            double b = coefficients[1];
            double c = coefficients[0];

            double discriminant = b * b - 4 * a * (c - binCenter_calibrated);

            if (discriminant < 0)
            {
                std::cerr << "Error: Negative discriminant. No real solution for bin center calibration." << std::endl;
                continue;
            }
            double sqrtDiscriminant = std::sqrt(discriminant);
            double x1 = (-b + sqrtDiscriminant) / (2 * a);
            double x2 = (-b - sqrtDiscriminant) / (2 * a);

            binCenter_original = (x1 >= 0) ? x1 : x2;
        }
        else
        {
            std::cerr << "Error: Unsupported polynomial degree. Supported degrees are 1 or 2." << std::endl;
            continue;
        }

        int bin_original = mainHist->GetXaxis()->FindBin(binCenter_original);

        if (bin_original < 1 || bin_original >= numBinsOriginal)
        {
            continue;
        }

        double interpolated_content = getInterpolatedContent(bin_original, binCenter_original);
        calibratedHist->SetBinContent(bin_calibrated, interpolated_content);
    }
}

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

void Histogram::setTotalArea()
{
    totalArea = 0;
    for (int bin = 1; bin <= mainHist->GetNbinsX(); ++bin)
    {
        totalArea += mainHist->GetBinContent(bin) * mainHist->GetBinWidth(bin);
    }
    std::cout << "Total area: " << totalArea << std::endl;
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
    std::cout << "Total area error: " << totalAreaError << std::endl;
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
    float ptError = pt * std::sqrt(
                             std::pow(areaPeakError / areaPeak, 2) + std::pow(totalAreaError / totalArea, 2));

    return ptError;
}

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

    for (int currentDegree = 1; currentDegree <= 2; ++currentDegree)
    {
        Eigen::MatrixXd X(n, currentDegree + 1);
        Eigen::VectorXd Y(n);

        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j <= currentDegree; ++j)
            {
                X(i, j) = std::pow(positions[i], j);
            }
            Y(i) = energies[i];
        }

        Eigen::VectorXd coeffs = X.colPivHouseholderQr().solve(Y);
        std::cout << "Polynomial degree: " << currentDegree << std::endl;
        std::cout << "Coefficients: " << coeffs[currentDegree] << std::endl;
        std::cout << "Threshold: " << polynomialFitThreshold << std::endl;
        if (std::abs(coeffs[currentDegree]) < polynomialFitThreshold && currentDegree > 1)
        {
            break;
        }

        degree = currentDegree;
        coefficients.clear();
        for (int i = 0; i <= currentDegree; ++i)
        {
            coefficients.push_back(coeffs(i));
        }
    }

    std::cout << "Final polynomial degree: " << degree << std::endl;
}
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
*/