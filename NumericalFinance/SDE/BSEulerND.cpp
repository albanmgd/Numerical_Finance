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

void BSEulerND::Simulate(double startTime, double endTime, size_t nbSteps, bool antitheticRV)
{
    SimBrownianND -> Simulate(startTime, endTime, nbSteps);
    double dt = (endTime - startTime) / nbSteps;

    /// we delete all existing paths
    for (SinglePath* path : Paths)
        delete path;
    Paths.clear();

    for (size_t d=0; d < Dimension; d++) {
        SinglePath *AssetPath = new SinglePath(startTime, endTime, nbSteps);
        double Spot = Spots[d];
        double Vol = Vols[d];
        AssetPath->AddValue(Spot);
        double lastInserted = Spot;

        //antithetic
        SinglePath *AntitheticPath = nullptr;
        double lastInsertedAntithetic = 0.0;
        if (antitheticRV){
            AntitheticPath = new SinglePath(startTime, endTime, nbSteps);
            AntitheticPath->AddValue(Spot);
            lastInsertedAntithetic = Spot;
        }

        for (size_t i=0; i < nbSteps; ++i) {
            double CorrelatedBM = SimBrownianND -> GetPath(d)->GetValue(i * dt);
            double nextValue = lastInserted
                               + lastInserted * (Rate * dt + Vol * CorrelatedBM);
            AssetPath->AddValue(nextValue);
            lastInserted = nextValue;

            //antithetic 
            if (antitheticRV){
                // we use -CorrelatedBM
                double nextValueAntithetic = lastInsertedAntithetic 
                + lastInsertedAntithetic * (Rate * dt - Vol * CorrelatedBM);
                AntitheticPath->AddValue(nextValueAntithetic);
                lastInsertedAntithetic = nextValueAntithetic;
            }
        }
        Paths.push_back(AssetPath);
        if (antitheticRV && AntitheticPath)
            Paths.push_back(AntitheticPath);
    }
}