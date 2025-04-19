#pragma once
#include "BlackScholesND.h"
#include "BrownianND.h"

class BSEulerND : public BlackScholesND
{
protected:
    int Dimension;
    vector<vector<double>>* Correls;
    BrownianND* SimBrownianND;
    bool HasToBeResimulated;
    vector<SinglePath*> AntitheticPaths;

public:
    BSEulerND(Normal* Gen, int dim, vector<double> spots, double rate, vector<double> vols, vector<vector<double>>* correls);
    void Simulate(double startTime, double endTime, size_t nbSteps, bool antitheticRV=false);
    ~BSEulerND();
};

