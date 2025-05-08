#pragma once

#include <iostream>
#include <vector>
#include "../RandomGenerator/Normal.h"

class BasketOption
{
protected:
    size_t Dimension;
    double K;
    double T;
    double Rate;
    std::vector<std::vector<double>>* Spots;
    std::vector<std::vector<double>>* Vols;
    std::vector<std::vector<double>>* Weights;
    std::vector<std::vector<double>>* Correls;


public:
    BasketOption(size_t dim, double K, double T, double Rate, std::vector<std::vector<double>>* Spots,
                 std::vector<std::vector<double>>* Vols, std::vector<std::vector<double>>* Weights,
                 std::vector<std::vector<double>>* Correls);
    virtual void PriceCall() = 0;
    ~BasketOption();
};

