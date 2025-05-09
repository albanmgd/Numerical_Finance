#pragma once
#include "BasketOption.h"

class EuropeanBasketOption : public BasketOption
{
protected:

public:
    EuropeanBasketOption(size_t dim, double K, double T, double Rate, std::vector<double> Spots,
                          std::vector<double> Vols, std::vector<double> Weights,
                          std::vector<std::vector<double>> Correls, Normal* Gen);
    ~EuropeanBasketOption();
    void PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic, bool UseControlVariate);
};
