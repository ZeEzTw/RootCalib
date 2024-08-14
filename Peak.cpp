#include "Peak.h"
#include <TF1.h>
#include <TH1D.h>
#include <cmath>
#include <iostream>

using namespace std;

Peak::Peak(TF1 *gausPeak, TH1D *hist)
    : position(0), gaus(nullptr), amplitude(0), sigma(0), area(0), leftLimit(0), rightLimit(0)
{
    if (gausPeak)
    {
        gaus = new TF1(*gausPeak); // Allocate and copy the function
        if (gaus)
        {
            position = gaus->GetParameter(1);
            amplitude = gaus->GetParameter(0);
            sigma = gaus->GetParameter(2);
            leftLimit = position - 3 * sigma;
            rightLimit = position + 3 * sigma;
            areaPeak(hist); // Assuming this is a member function that calculates the area
        }
    }
}

Peak::Peak(const Peak &other)
    : position(other.position), amplitude(other.amplitude), sigma(other.sigma), area(other.area),
      leftLimit(other.leftLimit), rightLimit(other.rightLimit)
{
    gaus = other.gaus ? new TF1(*other.gaus) : nullptr;
}

Peak &Peak::operator=(const Peak &other)
{
    if (this != &other)
    {
        position = other.position;
        amplitude = other.amplitude;
        sigma = other.sigma;
        area = other.area;
        leftLimit = other.leftLimit;
        rightLimit = other.rightLimit;

        delete gaus;
        gaus = other.gaus ? new TF1(*other.gaus) : nullptr;
    }
    return *this;
}

Peak::Peak(Peak &&other) noexcept
    : position(other.position), gaus(other.gaus), amplitude(other.amplitude), sigma(other.sigma),
      area(other.area), leftLimit(other.leftLimit), rightLimit(other.rightLimit)
{
    other.gaus = nullptr;
}

Peak &Peak::operator=(Peak &&other) noexcept
{
    if (this != &other)
    {
        position = other.position;
        amplitude = other.amplitude;
        sigma = other.sigma;
        area = other.area;
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
void Peak::areaPeak(TH1D *hist)
{
    double leftHeight = hist->GetBinContent(hist->FindBin(leftLimit - 3));
    double rightHeight = hist->GetBinContent(hist->FindBin(rightLimit + 3));
    if (leftHeight < 0.9 * rightHeight)
    {
        leftHeight = rightHeight;
    }
    else if (rightHeight < 0.9 * leftHeight)
    {
        rightHeight = leftHeight;
    }

    double backgroundHeight = (leftHeight + rightHeight) / 2;
    double width = abs(rightLimit - leftLimit);

    // Calculează aria sub gausiană în intervalul [leftLimit, rightLimit]
    double peakArea = gaus->Integral(leftLimit, rightLimit);
    double backgroundArea = backgroundHeight * width;

    // Calcularea ariei vârfului eliminând contribuția fundalului
    area = peakArea - backgroundArea;
}

float Peak::calculateResolution() const
{
    float FWHM = 2.3548 * sigma;
    float resolution = FWHM / amplitude;
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

    // Definirea limitelor MinDistance și MaxDistance
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
