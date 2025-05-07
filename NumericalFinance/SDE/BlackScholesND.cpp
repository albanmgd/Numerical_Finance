#include <string>
#include <stdexcept>
#include "BlackScholesND.h"
#include "../Utils/basic_functions.h"

BlackScholesND::BlackScholesND(Normal* Gen, vector<double> spots, double rate, vector<double> vols) :
        RandomProcess(Gen, 1), Spots(spots), Rate(rate), Vols(vols)
{
    size_t dim_spots = spots.size();
    size_t dim_vols = spots.size();
    if ((dim_spots != dim_vols) or (dim_spots == 0)){
        std::runtime_error("Please adjust the dimension of the inputs");
    }
    else
        Dimension = dim_spots;
}

BlackScholesND::~BlackScholesND(){}

double BlackScholesND::PriceBasketCallOption(double K, vector<double> Weights, double TimeToExpiry, vector<vector<double>> CorrelMatrix) {
    if (Weights.size() != Dimension || CorrelMatrix.size() != Dimension || CorrelMatrix[0].size() != Dimension)
        throw runtime_error("Dimension mismatch in weights or correlation matrix.");

    double X = 1.0;
    double adjDrift = 0.0;
    double sigmaEffSq = 0.0;

    // Compute X and adjDrift in one loop
    for (size_t i = 0; i < Dimension; ++i) {
        double w_i = Weights[i];
        double s_i = Spots[i];
        double v_i = Vols[i];

        X *= pow(s_i, w_i);
        adjDrift += w_i * v_i * v_i;

        for (size_t j = 0; j < Dimension; ++j) {
            double w_j = Weights[j];
            double v_j = Vols[j];
            sigmaEffSq += w_i * w_j * v_i * v_j * CorrelMatrix[i][j];
        }
    }

    double sigmaEff = sqrt(sigmaEffSq);
    double rEff = Rate - 0.5 * adjDrift + 0.5 * sigmaEffSq;

    double sqrtT = sqrt(TimeToExpiry);
    double d1 = (log(X / K) + (rEff + 0.5 * sigmaEffSq) * TimeToExpiry) / (sigmaEff * sqrtT);
    double d2 = d1 - sigmaEff * sqrtT;

    double callPrice = exp(-Rate * TimeToExpiry) *
                       (X * exp(rEff * TimeToExpiry) * normalCDF(d1) - K * normalCDF(d2));

    return callPrice;
}