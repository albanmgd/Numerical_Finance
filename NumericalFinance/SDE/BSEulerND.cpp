#include "BSEulerND.h"
#include "BrownianND.h"

BSEulerND::BSEulerND(Normal* Gen, int dim, vector <double> spots, double rate, vector<double> vols, vector<vector<double>>* correls) :
        BlackScholesND(Gen, spots, rate, vols), Dimension(dim), Correls(correls)
{
    SimBrownianND = new BrownianND(Gen, dim, correls);
    HasToBeResimulated = true;
}

BSEulerND::~BSEulerND()
{
}

void BSEulerND::Simulate(double startTime, double endTime, size_t nbSteps, bool antitheticRV)
{
    /*Main idea is: when we simulate the 'classic' paths, we also simulate antithetic paths if antitheticRV.
     * The next time we call this function, if antitheticRV then we return the previously simulated antithetic paths.
     * */

    // we delete all existing paths in every case since we'll either re-simulate them or re-write them
    for (SinglePath *path: Paths)
        delete path;
    Paths.clear();

    if (HasToBeResimulated) {
        SimBrownianND->Simulate(startTime, endTime, nbSteps);
        double dt = (endTime - startTime) / nbSteps;

        AntitheticPaths.clear();
        for (size_t d = 0; d < Dimension; d++) {
            SinglePath *AssetPath = new SinglePath(startTime, endTime, nbSteps);
            double Spot = Spots[d];
            double Vol = Vols[d];
            AssetPath->AddValue(Spot);
            double lastInserted = Spot;

            //antithetic
            SinglePath *AntitheticPath = nullptr;
            double lastInsertedAntithetic = 0.0;
            if (antitheticRV) {
                AntitheticPath = new SinglePath(startTime, endTime, nbSteps);
                AntitheticPath->AddValue(Spot);
                lastInsertedAntithetic = Spot;
            }

            for (size_t i = 0; i < nbSteps; ++i) {
                double CorrelatedBM = SimBrownianND->GetPath(d)->GetValue(i * dt);
                double nextValue = lastInserted
                                   + lastInserted * (Rate * dt + Vol * CorrelatedBM);
                AssetPath->AddValue(nextValue);
                lastInserted = nextValue;

                //antithetic
                if (antitheticRV) {
                    // we use -CorrelatedBM
                    double nextValueAntithetic = lastInsertedAntithetic
                                                 + lastInsertedAntithetic * (Rate * dt - Vol * CorrelatedBM);
                    AntitheticPath->AddValue(nextValueAntithetic);
                    lastInsertedAntithetic = nextValueAntithetic;
                }
            }
            Paths.push_back(AssetPath);
            if (antitheticRV) {
                AntitheticPaths.push_back(AntitheticPath);
            }
        }
        if (antitheticRV)
            HasToBeResimulated = false; /* could be in the loop but would lead to multiple affectations */
    }
    else if ((~HasToBeResimulated) & (antitheticRV)){
        for (SinglePath *AntitheticPath: AntitheticPaths){
            Paths.push_back(AntitheticPath);
        }
        HasToBeResimulated = true;
    }
}