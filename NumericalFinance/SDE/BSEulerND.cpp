#include "BSEulerND.h"
#include "BrownianND.h"

BSEulerND::BSEulerND(Normal* Gen, int dim, vector <double> spots, double rate, vector<double> vols, vector<vector<double>>* correls) :
        BlackScholesND(Gen, spots, rate, vols), Dimension(dim), Correls(correls)
{
    SimBrownianND = new BrownianND(Gen, dim, correls);
}

BSEulerND::~BSEulerND()
{
}

void BSEulerND::Simulate(double startTime, double endTime, size_t nbSteps)
{

    SimBrownianND -> Simulate(startTime, endTime, nbSteps);
    double dt = (endTime - startTime) / nbSteps;

    if(Paths.size() > 0) Paths.clear(); /* Deleting already existing paths if any */

    for (size_t d=0; d < Dimension; d++) {
        SinglePath *AssetPath = new SinglePath(startTime, endTime, nbSteps);
        double Spot = Spots[d];
        double Vol = Vols[d];
        AssetPath->AddValue(Spot);
        double lastInserted = Spot;
        for (size_t i=0; i < nbSteps; ++i) {
            double CorrelatedBM = SimBrownianND -> GetPath(d)->GetValue(i * dt);
            double nextValue = lastInserted
                               + lastInserted * (Rate * dt + Vol * CorrelatedBM);
            AssetPath->AddValue(nextValue);
            lastInserted = nextValue;
        }
        Paths.push_back(AssetPath);
    }
}