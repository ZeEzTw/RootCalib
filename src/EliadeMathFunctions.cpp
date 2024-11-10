#include "../include/EliadeMathFunctions.h"
#include <cmath>
#include <algorithm>

std::vector<std::vector<double>> EliadeMathFunctions::multiplyTransposeMatrix(const std::vector<std::vector<double>> &X) {
    int rows = X.size();
    int cols = X[0].size();
    std::vector<std::vector<double>> XtX(cols, std::vector<double>(cols, 0.0));

    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < rows; ++k) {
                sum += X[k][i] * X[k][j];
            }
            XtX[i][j] = XtX[j][i] = sum;
        }
    }

    return XtX;
}

std::vector<double> EliadeMathFunctions::multiplyTransposeVector(const std::vector<std::vector<double>> &X, const std::vector<double> &Y) {
    int rows = X.size();
    int cols = X[0].size();
    std::vector<double> XtY(cols, 0.0);

    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            XtY[i] += X[j][i] * Y[j];
        }
    }

    return XtY;
}

std::vector<double> EliadeMathFunctions::solveSystem(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
    int n = A.size();
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(n + 1));

    // Augment the matrix A with the vector b
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = b[i];
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        double maxElement = std::abs(augmentedMatrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > maxElement) {
                maxElement = std::abs(augmentedMatrix[k][i]);
                maxRow = k;
            }
        }

        if (maxRow != i) {
            std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);
        }

        for (int k = i + 1; k < n; ++k) {
            double c = -augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j <= n; ++j) {
                if (i == j) {
                    augmentedMatrix[k][j] = 0;
                } else {
                    augmentedMatrix[k][j] += c * augmentedMatrix[i][j];
                }
            }
        }
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; --k) {
            augmentedMatrix[k][n] -= augmentedMatrix[k][i] * x[i];
        }
    }

    return x;
}
