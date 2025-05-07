#include "basic_functions.h"
#include <cmath>

double normalCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
};
