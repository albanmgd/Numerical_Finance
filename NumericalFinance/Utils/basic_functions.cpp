#include "basic_functions.h"
#include <cmath>
#include <iostream>
#include <vector>

double normalCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
};

double meanVector(std::vector<double> vec){
    double nbElements = vec.size();
    double sum = 0.0;
    for (size_t i=0; i<nbElements; i++){
        sum += vec[i];
    }
    return sum / nbElements;
}

double varianceVector(std::vector<double> vec){
    double nbElements = vec.size();
    double mean = meanVector(vec);
    double var = 0.0;
    for (size_t i=0; i<nbElements; i++){
        var += pow(vec[i] - mean, 2);
    }
    return var / nbElements;
}

double laguerrePolynomial(std::size_t l, double x) {
    if (l == 0) return 1.0;
    if (l == 1) return 1.0 - x;

    double L_prev2 = 1.0;          // L_0(x)
    double L_prev1 = 1.0 - x;      // L_1(x)
    double L_curr;

    for (size_t n = 2; n <= l; ++n) {
        L_curr = ((2 * n - 1 - x) * L_prev1 - (n - 1) * L_prev2) / n;
        L_prev2 = L_prev1;
        L_prev1 = L_curr;
    }

    return L_curr;
}