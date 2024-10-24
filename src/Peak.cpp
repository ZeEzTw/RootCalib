#include "../include/Peak.h"
#include <TF1.h>
#include <TH1D.h>
#include <cmath>
#include <iostream>

using namespace std;

Peak::Peak(TF1 *gausPeak, TH1D *hist)
    : position(0), associatedPosition(0), gaus(nullptr), amplitude(0), sigma(0), area(0), areaError(0), leftLimit(0), rightLimit(0)
{
    if (gausPeak)
    {
        gaus = new TF1(*gausPeak);
        if (gaus)
        {
            position = gaus->GetParameter(1);
            amplitude = gaus->GetParameter(0);
            sigma = gaus->GetParameter(2);
            leftLimit = position - 3 * sigma;
            rightLimit = position + 3 * sigma;
            areaPeak(hist);
        }
    }
}

Peak::Peak(const Peak &other)
    : position(other.position), associatedPosition(other.associatedPosition), amplitude(other.amplitude), sigma(other.sigma),
      area(other.area), areaError(other.areaError), leftLimit(other.leftLimit), rightLimit(other.rightLimit)
{
    gaus = other.gaus ? new TF1(*other.gaus) : nullptr;
}

Peak &Peak::operator=(const Peak &other)
{
    if (this != &other)
    {
        position = other.position;
        amplitude = other.amplitude;
        associatedPosition = other.associatedPosition;
        sigma = other.sigma;
        area = other.area;
        areaError = other.areaError;
        leftLimit = other.leftLimit;
        rightLimit = other.rightLimit;

        delete gaus;
        gaus = other.gaus ? new TF1(*other.gaus) : nullptr;
    }
    return *this;
}

Peak::Peak(Peak &&other) noexcept
    : position(other.position), associatedPosition(other.associatedPosition), gaus(other.gaus), amplitude(other.amplitude), sigma(other.sigma),
      area(other.area), areaError(other.areaError), leftLimit(other.leftLimit), rightLimit(other.rightLimit)
{
    other.gaus = nullptr;
}

Peak &Peak::operator=(Peak &&other) noexcept
{
    if (this != &other)
    {
        position = other.position;
        associatedPosition = other.associatedPosition;
        amplitude = other.amplitude;
        sigma = other.sigma;
        area = other.area;
        areaError = other.areaError;
        leftLimit = other.leftLimit;
        rightLimit = other.rightLimit;

        delete gaus;
        gaus = other.gaus;
        other.gaus = nullptr;
    }
    return *this;
}

Peak::~Peak()
{
    delete gaus;
}
void Peak::setPosition(double pos)
{
    position = pos;
}

double Peak::getPosition() const
{
    return position;
}

void Peak::setAssociatedPosition(double pos)
{
    associatedPosition = pos;
}

double Peak::getAssociatedPosition() const
{
    return associatedPosition;
}

void Peak::setAmplitude(double amp)
{
    amplitude = amp;
}

double Peak::getAmplitude() const
{
    return amplitude;
}

void Peak::setSigma(double sig)
{
    sigma = sig;
}

double Peak::getSigma() const
{
    return sigma;
}
double Peak::getFWHM() const
{
    return 2.3548 * sigma;
}
double Peak::getMean() const
{
    return gaus->GetParameter(1);
}
void Peak::createGaussianFunction()
{
    // Funcția poate fi implementată dacă este necesar
}

TF1 *Peak::getGaussianFunction() const
{
    return gaus;
}
void Peak::setArea(double a)
{
    area = a;
}
double Peak::getArea() const
{
    return area;
}

double Peak::getAreaError() const
{
    return areaError;
}

void Peak::setLeftLimit(float left)
{
    leftLimit = left;
}
float Peak::getLeftLimit() const
{
    return leftLimit;
}
void Peak::setRightLimit(float right)
{
    rightLimit = right;
}
float Peak::getRightLimit() const
{
    return rightLimit;
}
/*void Peak::areaPeak(TH1D *hist)
{
    // Assume backgroundModel has been fitted elsewhere and provides background at each bin
    TF1 *backgroundModel = ...; // This should be defined based on your background fitting

    double peakArea = 0.0;

    // Loop over bins within the peak region
    for (int bin = hist->FindBin(leftLimit); bin <= hist->FindBin(rightLimit); ++bin) {
        double binContent = hist->GetBinContent(bin);
        double binCenter = hist->GetBinCenter(bin);
        double background = backgroundModel->Eval(binCenter);

        peakArea += (binContent - background) * hist->GetBinWidth(bin);
    }

    area = peakArea;
}*/

void Peak::areaPeak(TH1D *hist)
{
    int leftBin = hist->FindBin(leftLimit);
    int rightBin = hist->FindBin(rightLimit);

    double leftHeight = hist->GetBinContent(leftBin);
    double rightHeight = hist->GetBinContent(rightBin);

    double errorLeft = hist->GetBinError(leftBin);
    double errorRight = hist->GetBinError(rightBin);

    double peakArea = 0.0;
    double peakAreaErrorSq = 0.0;

    double backgroundArea = 0.0;
    double backgroundAreaErrorSq = 0.0;

    for (int bin = leftBin; bin <= rightBin; ++bin)
    {
        double binCenter = hist->GetBinCenter(bin);
        double binWidth = hist->GetBinWidth(bin);
        double binContent = hist->GetBinContent(bin);
        double binError = hist->GetBinError(bin);

        double fraction = (binCenter - leftLimit) / (rightLimit - leftLimit);
        double background = leftHeight + (rightHeight - leftHeight) * fraction;

        if (binContent <= 0 || binError < 0)
        {
            //std::cerr << "Invalid bin content or bin error at bin " << bin << std::endl;
            continue;
        }

        double backgroundError = std::sqrt(std::pow(errorLeft, 2) + std::pow(errorRight, 2)) * fraction;

        peakArea += (binContent * binWidth);
        peakAreaErrorSq += std::pow(binError * binWidth, 2);

        backgroundArea += (background * binWidth);
        backgroundAreaErrorSq += std::pow(backgroundError * binWidth, 2);
    }

    area = peakArea - backgroundArea;

    double totalErrorSq = abs(peakAreaErrorSq + backgroundAreaErrorSq);
    areaError = std::sqrt(totalErrorSq);
}

double Peak::calculateResolutionError() const
{
    // Obținem erorile pentru amplitudine și sigma
    double amplitudeError = gaus->GetParError(0); // Eroarea parametrului de amplitudine
    double sigmaError = gaus->GetParError(2);     // Eroarea parametrului de sigma

    // Calculăm FWHM
    double FWHM = 2.3548 * sigma;

    // Calculăm derivatele pentru rezoluție
    double dR_dA = -FWHM / (amplitude * amplitude); // Derivata rezoluției în raport cu amplitudinea
    double dR_dSigma = 2.3548 / amplitude;          // Derivata rezoluției în raport cu sigma

    // Calculăm eroarea rezoluției folosind propagarea erorilor
    double resolutionError = std::sqrt(std::pow(dR_dA * amplitudeError, 2) +
                                       std::pow(dR_dSigma * sigmaError, 2));

    // Returnăm eroarea rezoluției
    return resolutionError;
}

double Peak::calculateResolution() const
{
    double FWHM = 2.3548 * sigma;
    double resolution = FWHM / amplitude;
    return resolution;
}

void Peak::findStartOfPeak(TH1D *hist, int maxBin, double &leftLimitPosition, double &rightLimitPosition)
{
    double left = abs(leftLimitPosition - maxBin);
    double right = abs(rightLimitPosition - maxBin);

    if (left > right)
    {
        right = left;
    }
    else
    {
        left = right;
    }

    double MinDistance = 1.0;  // Exemplar
    double MaxDistance = 10.0; // Exemplar

    if (left > MaxDistance || left < MinDistance)
    {
        leftLimitPosition = maxBin - MinDistance;
    }
    if (right > MaxDistance || right < MinDistance)
    {
        rightLimitPosition = maxBin + MinDistance;
    }
}
void Peak::outputDataJson(std::ofstream &file) const
{
    file << "{\n";
    file << "\t\"position\": " << position << ",\n";
    file << "\t\"amplitude\": " << amplitude << ",\n";
    file << "\t\"sigma\": " << sigma << ",\n";
    file << "\t\"area\": " << area << ",\n";
    file << "\t\"leftLimit\": " << leftLimit << ",\n";
    file << "\t\"rightLimit\": " << rightLimit << "\n";
    file << "}";
}
