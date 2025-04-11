#pragma once
#include "RandomProcess.h"

class BlackScholesND : public RandomProcess
{
protected:
    vector<double> Spots;
    double Rate;
    vector<double> Vols;

public:
    BlackScholesND(Normal* Gen, vector<double> spot, double rate, vector<double> vol);
    ~BlackScholesND();
};
