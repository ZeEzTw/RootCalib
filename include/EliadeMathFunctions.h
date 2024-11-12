/**
 * @brief Specialized mathematical utility class for histogram processing
 * 
 * This class provides essential mathematical operations required for histogram analysis
 * and data processing, focusing on matrix operations commonly used in spectroscopy
 * and statistical analysis of histogram data.
 * 
 * Core functionalities:
 * - Matrix operations for histogram data transformation
 * - Transpose multiplication for data correlation
 * - Linear system solving for calibration calculations
 * 
 * Extensibility notes:
 * - Functions can be modified from static to member functions for an object-oriented approach
 * - Additional mathematical operations can be added as needed
 * - The current implementation uses double, but can be adapted using templates for other data types
 * 
 * Example usage in histogram processing:
 *     std::vector<std::vector<double>> histogramData = ...;
 *     auto correlationMatrix = EliadeMathFunctions::multiplyTransposeMatrix(histogramData);
 */

#ifndef ELIADEMATHFUNCTIONS_H
#define ELIADEMATHFUNCTIONS_H

#include <vector>

class EliadeMathFunctions {
public:
    static std::vector<std::vector<double>> multiplyTransposeMatrix(const std::vector<std::vector<double>> &X);
    static std::vector<double> multiplyTransposeVector(const std::vector<std::vector<double>> &X, const std::vector<double> &Y);
    /**
     * @brief Solves a system of linear equations.
     * 
     * This function takes a matrix A and a vector b, and solves the system of linear equations A * x = b.
     * 
     * @param A The coefficient matrix of the system of equations.
     * @param b The right-hand side vector of the system of equations.
     * @return A vector containing the solution to the system of equations.
     */
    static std::vector<double> solveSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b);
};

#endif // ELIADEMATHFUNCTIONS_H
