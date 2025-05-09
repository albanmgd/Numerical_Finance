#pragma once

#include <iostream>
#include <vector>
#include "../RandomGenerator/Normal.h"

class BasketOption
{
protected:
    size_t Dimension; // nb of assets in the basket
    double K; // strike
    double T; // expiry
    double Rate;
    std::vector<double> Spots;
    std::vector<double> Vols;
    std::vector<double> Weights; // weights associated to the stocks in the basket
    std::vector<std::vector<double>> Correls; // Correlation matrix
    Normal* Generator;

public:
    BasketOption(size_t dim, double K, double T, double Rate, std::vector<double> Spots,
                 std::vector<double> Vols, std::vector<double> Weights,
                 std::vector<std::vector<double>> Correls, Normal* Gen);
    virtual std::vector<double> PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic,
                                          bool UseControlVariate) = 0;
    ~BasketOption();
};

