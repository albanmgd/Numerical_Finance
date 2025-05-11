#pragma once
#include "BasketOption.h"

class BermudanBasketOption : public BasketOption
{
protected:
    size_t L;

public:
    BermudanBasketOption(size_t Dim, double K, double T, double Rate, std::vector<double> Spots,
                         std::vector<double> Vols, std::vector<double> Weights,
                         std::vector<std::vector<double>> Correls, Normal* Gen, size_t L);
    ~BermudanBasketOption();

    std::vector<double> PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic,
                                  bool UseControlVariate);

};
