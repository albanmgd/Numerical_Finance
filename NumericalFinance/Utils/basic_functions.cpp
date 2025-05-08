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
