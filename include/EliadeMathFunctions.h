#ifndef ELIADEMATHFUNCTIONS_H
#define ELIADEMATHFUNCTIONS_H

#include <vector>

class EliadeMathFunctions {
public:
    // Multiplies a matrix by its transpose
    static std::vector<std::vector<double>> multiplyTransposeMatrix(const std::vector<std::vector<double>> &X);

    // Multiplies the transpose of a matrix by a vector
    static std::vector<double> multiplyTransposeVector(const std::vector<std::vector<double>> &X, const std::vector<double> &Y);

    // Solves a linear system using Gaussian elimination
    static std::vector<double> solveSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b);
};

#endif // ELIADEMATHFUNCTIONS_H
