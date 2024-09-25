#ifndef PEAK_H
#define PEAK_H

#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>
class Peak
{
private:
    double position;
    TF1 *gaus;
    double amplitude;
    double sigma;
    double area;
    float leftLimit, rightLimit;

public:
    Peak(TF1 *gausPeak, TH1D *hist); // Constructor
    // Constructor with TF1 parameter
    Peak(TF1 *gausPeak);

    // Copy constructor
    Peak(const Peak &other);

    // Copy assignment operator
    Peak &operator=(const Peak &other);

    // Move constructor
    Peak(Peak &&other) noexcept;

    // Move assignment operator
    Peak &operator=(Peak &&other) noexcept;

    // Destructor
    ~Peak();

    void setPosition(double pos);
    double getPosition() const;

    void setAmplitude(double amp);
    double getAmplitude() const;

    void setSigma(double sig);
    double getSigma() const;

    void setArea(double a);
    double getArea() const;

    void setLeftLimit(float left);
    float getLeftLimit() const;

    void setRightLimit(float right);
    float getRightLimit() const;

    void outputDataJson(std::ofstream &file) const;
    void createGaussianFunction();
    TF1 *getGaussianFunction() const;
    double getFWHM() const;
    double getMean() const;
    void areaPeak(TH1D *hist);
    float calculateResolution() const;
    void findStartOfPeak(TH1D *hist, int maxBin, double &leftLimitPosition, double &rightLimitPosition);
};

#endif // PEAK_H
       // sa fac functia sa se poata redefinii un peak, adica sa mearga sa zica peak 6 e de fapt la pozitia asta, sa mearga acolo si sa faca fitul