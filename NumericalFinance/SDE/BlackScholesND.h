#pragma once
#include "RandomProcess.h"

class BlackScholesND : public RandomProcess
{
protected:
    vector<double> Spots;
    double Rate;
    vector<double> Vols;
    size_t Dimension;

public:
    BlackScholesND(Normal* Gen, vector<double> spot, double rate, vector<double> vol);
    ~BlackScholesND();
    double PriceBasketCallOption(double K, vector<double> weights, double TimeToExpiry, vector<vector<double>> CorrelMatrix);
};
