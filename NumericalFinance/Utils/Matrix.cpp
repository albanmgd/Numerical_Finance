#include "Matrix.h"
#include <cmath>
#include <stdexcept>

// Default constructor - empty 2*2 matrix
Matrix::Matrix(): data(2, std::vector<double>(2, 0.0)) {}

// Parameterized constructor
Matrix::Matrix(const std::vector<std::vector<double>>& vect): data(vect) {}

void Matrix::setMatrixData(const Matrix& other_mat) {
    this->data = other_mat.data;
}

void Matrix::addRow(const std::vector<double>& newRow, int newRowIndex) {
    if (newRowIndex >= 0 && newRowIndex <= data.size()) { // Ensure n is within valid bounds
        data.insert(data.begin() + newRowIndex, newRow);
    } else {
        std::runtime_error("Index out of bounds!");
    }
}

Matrix Matrix::getTranspose() {
    // Creating a matrix of the same size
    std::vector<std::vector<double>> solution(this->data[0].size(), std::vector<double> (this->data.size()));

    // Filling solution-matrix
    for(size_t i = 0; i < this->data.size(); i++) {
        for(size_t j = 0; j < this->data[0].size(); j++) {
            solution[j][i] = this->data[i][j];
        }
    }
    return Matrix(solution);
}

bool Matrix::luDecomposition(const std::vector<std::vector<double>>& matrix,
                             std::vector<std::vector<double>>& L,
                             std::vector<std::vector<double>>& U) {
    size_t n = matrix.size();

    // Initialize L and U
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        // Upper triangular matrix U
        for (size_t k = i; k < n; ++k) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = matrix[i][k] - sum;
        }

        // Lower triangular matrix L
        for (size_t k = i; k < n; ++k) {
            if (i == k) {
                L[i][i] = 1.0; // Diagonal of L is 1
            } else {
                double sum = 0.0;
                for (size_t j = 0; j < i; ++j) {
                    sum += L[k][j] * U[j][i];
                }
                if (std::fabs(U[i][i]) < 1e-9) {
                    return false; // Singular matrix
                }
                L[k][i] = (matrix[k][i] - sum) / U[i][i];
            }
        }
    }

    return true;
}

double Matrix::determinantLU() {
    size_t n = data.size();
    std::vector<std::vector<double>> L, U;

    if (!luDecomposition(data, L, U)) {
        return 0.0; // Singular matrix
    }

    double det = 1.0;
    for (size_t i = 0; i < n; ++i) {
        det *= U[i][i];
    }

    return det;
}

Matrix Matrix::inverseLU() {
    size_t n = data.size();
    std::vector<std::vector<double>> L, U;

    if (!luDecomposition(data, L, U)) {
        throw std::runtime_error("Matrix is singular, inverse does not exist.");
    }
    // Initialize identity matrix
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        I[i][i] = 1.0;
    }
    // Forward substitution to solve L * Y = I
    std::vector<std::vector<double>> Y(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L[i][k] * Y[k][j];
            }
            Y[i][j] = I[i][j] - sum;
        }
    }
    // Back substitution to solve U * X = Y
    std::vector<std::vector<double>> X(n, std::vector<double>(n, 0.0));
    for (int i = n - 1; i >= 0; --i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = i + 1; k < n; ++k) {
                sum += U[i][k] * X[k][j];
            }
            X[i][j] = (Y[i][j] - sum) / U[i][i];
        }
    }
    Matrix InvMat = Matrix(X);
    return InvMat;
}

Matrix Matrix::operator-(const Matrix& other) const {
    // Check nb of rows and cols
    if (data.size() != other.data.size()) {
        throw std::runtime_error("Matrices don't have the same number of rows.");
    }
    if (data[0].size() != other.data[0].size()) {
        throw std::runtime_error("Matrices don't have the same number of columns.");
    }
    std::vector<std::vector<double>> diff_vect(data.size(), std::vector<double>(data[0].size()));
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            diff_vect[i][j] = data[i][j] - other.data[i][j];
        }
    }
    return Matrix(diff_vect);
}

Matrix Matrix::operator+(const Matrix& other) const {
    // Check nb of rows and cols
    if (data.size() != other.data.size()) {
        throw std::runtime_error("Matrices don't have the same number of rows.");
    }
    if (data[0].size() != other.data[0].size()) {
        throw std::runtime_error("Matrices don't have the same number of columns.");
    }
    std::vector<std::vector<double>> sum_vect(data.size(), std::vector<double>(data[0].size()));
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            sum_vect[i][j] = data[i][j] + other.data[i][j];
        }
    }
    return Matrix(sum_vect);
}

Matrix Matrix::operator*(const Matrix& other) const {
    // Ensure matrices are non-empty
    if (data.empty() || other.data.empty() || other.data[0].empty()) {
        throw std::runtime_error("Matrices cannot be multiplied: one or both matrices are empty.");
    }

    // Check dimension compatibility
    if (data[0].size() != other.data.size()) {
        throw std::runtime_error("Matrices cannot be multiplied: number of columns in the first matrix != number of rows in the second matrix.");
    }

    // Initialize the result matrix
    size_t rows = data.size();
    size_t cols = other.data[0].size();
    size_t common_dim = other.data.size();
    std::vector<std::vector<double>> prod_vect(rows, std::vector<double>(cols, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < common_dim; ++k) {
                prod_vect[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    return Matrix(prod_vect);
}