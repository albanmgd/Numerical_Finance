# pragma once
#include "RandomProcess.h"

class BrownianND: public RandomProcess
{
protected:
    int Dimension;
    vector < vector<double> >* CorrelationMatrix;
    vector < vector<double> > CholeskyDecomposition;

public:
    BrownianND (RandomGenerator* Gen , int dim, vector <vector<double>>* Corr);
    std::vector<std::vector<double>> getCholeskyDecomposition( vector < vector<double> >* Corr);
    void Simulate (double startTime, double endTime, size_t nbSteps, bool antitheticRV=false);
};
