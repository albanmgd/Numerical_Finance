#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

class PAdic {
private:
    int Base;
    int Precision;

public:
    PAdic(int p, int prec = 20);
    std::vector<int> double_to_padic(double x); /* Converts a double to its p-adic decomposition */
    double padic_to_double(std::vector<int>* digits); /* Converts a vector of int representing the p-adic decomp. to a double */
    double add(double firstNb, double secondNb); /* p-adic summation on two doubles */
//    std::vector<int> rotate(double* number, double* angle); /* rotation of a given angle */

};

