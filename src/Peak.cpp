#include "../include/Peak.h"
#include <TF1.h>
#include <TH1D.h>
#include <cmath>
#include <iostream>

using namespace std;

Peak::Peak(TF1 *gausPeak, TH1D *hist)
    : position(0), associatedPosition(0), gaus(nullptr), amplitude(0), sigma(0), area(0), areaError(0), leftLimit(0), rightLimit(0)
{
    if (!gausPeak) return;
    
    gaus = new TF1(*gausPeak);
    if (!gaus) return;
    
    position = gaus->GetParameter(1);
    amplitude = gaus->GetParameter(0);
    sigma = gaus->GetParameter(2);
    leftLimit = position - SIGMA_MULTIPLIER * sigma;
    rightLimit = position + SIGMA_MULTIPLIER * sigma;
    
    if (hist) areaPeak(hist);
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

void Peak::areaPeak(TH1D* hist) 
{
    if (!hist) return;
    
    int leftBin = hist->FindBin(leftLimit);
    int rightBin = hist->FindBin(rightLimit);
    
    double leftHeight = hist->GetBinContent(leftBin);
    double rightHeight = hist->GetBinContent(rightBin);
    double errorLeft = hist->GetBinError(leftBin);
    double errorRight = hist->GetBinError(rightBin);
    
    double totalArea = 0.0;
    double totalError = 0.0;
    double bgArea = 0.0;
    double bgError = 0.0;
    
    for (int bin = leftBin; bin <= rightBin; ++bin) {
        double binContent = hist->GetBinContent(bin);
        double binError = hist->GetBinError(bin);
        
        if (!isValidBin(binContent, binError)) continue;
        
        double binWidth = hist->GetBinWidth(bin);
        double binCenter = hist->GetBinCenter(bin);
        
        double fraction = (binCenter - leftLimit) / (rightLimit - leftLimit);
        double background = leftHeight + (rightHeight - leftHeight) * fraction;
        
        totalArea += binContent * binWidth;
        totalError += binError * binError * binWidth * binWidth;
        bgArea += background * binWidth;
        bgError += (errorLeft * errorLeft + errorRight * errorRight) * 
                   fraction * fraction * binWidth * binWidth;
    }
    
    area = totalArea - bgArea;
    areaError = std::sqrt(std::abs(totalError + bgError));
}

double Peak::calculateResolutionError() const
{
    if (!gaus) return 0.0;
    
    double amplitudeError = gaus->GetParError(0);
    double sigmaError = gaus->GetParError(2);
    double fwhm = getFWHM();
    
    double dR_dA = -fwhm / (amplitude * amplitude);
    double dR_dSigma = FWHM_CONSTANT / amplitude;
    
    return std::sqrt(dR_dA * dR_dA * amplitudeError * amplitudeError + 
                    dR_dSigma * dR_dSigma * sigmaError * sigmaError);
}

double Peak::calculateResolution() const 
{
    return getFWHM() / amplitude;
}

void Peak::findStartOfPeak(TH1D *hist, int maxBin, double &leftLimitPosition, double &rightLimitPosition)
{
    double left = std::abs(leftLimitPosition - static_cast<double>(maxBin));
    double right = std::abs(rightLimitPosition - static_cast<double>(maxBin));

    if (left > right)
    {
        right = left;
    }
    else
    {
        left = right;
    }

    double MinDistance = Constants::MIN_DISTANCE;
    double MaxDistance = Constants::MAX_DISTANCE;

    if (left > MaxDistance || left < MinDistance)
    {
        leftLimitPosition = static_cast<double>(maxBin) - MinDistance;
    }
    if (right > MaxDistance || right < MinDistance)
    {
        rightLimitPosition = static_cast<double>(maxBin) + MinDistance;
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
