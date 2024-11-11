/**
 * @class Peak
 * @brief A class for analyzing and storing properties of histogram peaks
 *
 * The Peak class manages the properties and analysis of peaks in histograms.
 * Key features:
 * - Stores peak properties (position, amplitude, sigma, area)
 * - Calculates peak area with background subtraction
 * - Computes resolution and FWHM (Full Width at Half Maximum)
 * - Handles Gaussian fitting functions
 * - Provides JSON output capabilities
 * - Supports copy/move semantics
 *
 * Note: Some methods are not used in the codebase (ex: findStartOfPeak()) but are retained to facilitate future
 * extensions, allowing the class to perform different calculations on peaks or to
 * handle additional information.
 */

#ifndef PEAK_H
#define PEAK_H

#include <TF1.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

class Peak {
private:
    // Constants
    enum Constants {
        SIGMA_MULTIPLIER = 3,
        MIN_DISTANCE = 1,
        MAX_DISTANCE = 10
    };
    static constexpr double FWHM_CONSTANT = 2.3548;

    // Core peak properties
    double position;
    double associatedPosition;
    double amplitude;
    double sigma;
    
    // Area related
    double area;
    double areaError;
    float leftLimit;
    float rightLimit;
    
    // Gaussian function
    TF1* gaus;

    // Helper methods
    bool isValidBin(double content, double error) const {
        return content > 0 && error >= 0;
    }

public:
    // Constructors and destructor
    explicit Peak(TF1* gausPeak, TH1D* hist = nullptr);
    Peak(const Peak& other);
    Peak& operator=(const Peak& other);
    Peak(Peak&& other) noexcept;
    Peak& operator=(Peak&& other) noexcept;
    ~Peak();

    // Analysis methods
    void areaPeak(TH1D* hist);
    double getFWHM() const { return FWHM_CONSTANT * sigma; }
    double getMean() const {return gaus->GetParameter(1); };
    double calculateResolution() const;
    double calculateResolutionError() const;
    //not used in code
    void findStartOfPeak(TH1D* hist, int maxBin, double& leftLimitPosition, double& rightLimitPosition);

    // Utility methods
    TF1* getGaussianFunction() const { return gaus; }
    void createGaussianFunction();
    void outputDataJson(std::ofstream& file) const;

    // Getters
    double getPosition() const { return position; }
    double getAssociatedPosition() const { return associatedPosition; }
    double getAmplitude() const { return amplitude; }
    double getSigma() const { return sigma; }
    double getArea() const { return area; }
    double getAreaError() const { return areaError; }
    float getLeftLimit() const { return leftLimit; }
    float getRightLimit() const { return rightLimit; }

    // Setters
    void setPosition(double pos) { position = pos; }
    void setAssociatedPosition(double pos) { associatedPosition = pos; }
    void setAmplitude(double amp) { amplitude = amp; }
    void setSigma(double sig) { sigma = sig; }
    void setArea(double a) { area = a; }
    void setLeftLimit(float left) { leftLimit = left; }
    void setRightLimit(float right) { rightLimit = right; }
};

#endif // PEAK_H