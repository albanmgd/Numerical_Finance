#include "BermudanBasketOption.h"
#include "../Utils/Matrix.h"
#include "../Utils/basic_functions.h"
#include "../SDE/BSEulerND.h"
#include <memory>

BermudanBasketOption::BermudanBasketOption(
        size_t Dim, double K, double T, double Rate, std::vector<double> Spots,
        std::vector<double> Vols, std::vector<double> Weights,
        std::vector<std::vector<double>> Correls, Normal* Gen, size_t L) :
        BasketOption(Dim, K, T, Rate, Spots, Vols,Weights,Correls, Gen), L(L)
        {};

std::vector<double> BermudanBasketOption::PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic,
                                                    bool UseControlVariate){
    cout << "Starting the MC Simulation for the Bermudean Call..." << endl;
    clock_t start, end;
    start = clock();

    BSEulerND TestScheme = BSEulerND(Generator, Dimension, Spots, Rate, Vols, &Correls);
    // Computing the theoretical price in case we use the control variate - avoids multiple if tests
    double TheoreticalPrice = TestScheme.PriceBasketCallOption(K, Weights, T, Correls);

    /* Begin by generating directly nbSim asset paths for the d assets & computing the basket values associated */
    double delta_t = T / NbSteps;
    vector<vector<std::unique_ptr<SinglePath>>> AllAssetPaths(NbSims);
    std::vector<SinglePath*> basketValues(NbSims, nullptr);

    for (size_t nSimul=0; nSimul < NbSims; nSimul++) {
        AllAssetPaths[nSimul].resize(Dimension); /* resizing to hold d SinglePath */
        TestScheme.Simulate(0, T, NbSteps, UseAntithetic); /* Simulating the paths - enables us to reuse antithetic control variate*/
        for (size_t nbAsset=0; nbAsset < Dimension; nbAsset++) {
            AllAssetPaths[nSimul][nbAsset] = std::make_unique<SinglePath>(*TestScheme.GetPath(nbAsset)); /* need to do copy of paths otherwise pointers get deleted every time we re-call the generate function*/
        }
        /* Can finally compute the basket values for that simulation */
        SinglePath* basketValueSimulation = new SinglePath(0.0, T, NbSteps);
        for (size_t i=0; i<NbSteps; i++){
            double currentTimestep = i * delta_t;
            double basketValueSim = 0.0;
            for (size_t d=0; d<Dimension; d++){
                basketValueSim += Weights[d] * AllAssetPaths[nSimul][d]->GetValue(currentTimestep);
            }
            basketValueSimulation->AddValue(basketValueSim);
        }
        basketValues[nSimul] = basketValueSimulation;
    }
    // Backward induction for the stopping times
    std::vector<std::vector<double>> stoppingTimes(NbSims, std::vector<double>(NbSteps, 0.0));
    for (size_t n = 0; n < NbSims; ++n) {
        stoppingTimes[n][NbSteps - 1] = 1.0;
    }

    // Looping through time
    for (size_t i=(NbSteps - 2); i >= 1; i--){
        std::vector<std::vector<double>> basisVectors(NbSims, std::vector<double>(L, 0.0));
        double currentTimestep = i * delta_t; /* going backward */
        // Initializing the right hand part of the sum to optimize - colum vector
        vector<vector<double>> basisDecompositionVector (L, vector<double>(1, 0.0));

        /*vector<double> basketValuesAtTimestep(NbSims, 0.0);
        for (size_t nSimul=0; nSimul < NbSims; nSimul++) {
            basketValuesAtTimestep[nSimul] = basketValues[nSimul] ->GetValue(currentTimestep);
        }
        double crossSectionalMeanBasketValue = meanVector(basketValuesAtTimestep);
        double crossSectionalStdBasketValue = sqrt(varianceVector(basketValuesAtTimestep));*/

        // Going through all the simulations for one time-step
        for (size_t nSimul=0; nSimul < NbSims; nSimul++) {
            double previousStoppingTime = stoppingTimes[nSimul][i + 1];
            double currentBasketValue = basketValues[nSimul]->GetValue(currentTimestep);
            double pastStoppingTimeBasketValue = basketValues[nSimul]->GetValue(previousStoppingTime);
            double mulFactor = exp(- Rate * (previousStoppingTime - currentTimestep)) * max<double>(pastStoppingTimeBasketValue - K, 0.0);
            /* Constructing the P vector */
            for (size_t l=0; l < L; l++){
                //double basisScalar = laguerrePolynomial(l, (currentBasketValue - crossSectionalMeanBasketValue)/crossSectionalStdBasketValue);
                double basisScalar = laguerrePolynomial(l, currentBasketValue);
                basisVectors[nSimul][l] = basisScalar;
                basisDecompositionVector[l][0] += mulFactor * basisScalar;
            }
        }
        // Can now construct the Phi Matrix
        Matrix Phi = Matrix(basisVectors);
        Matrix H = Matrix(Phi.getTranspose() * Phi);

        Matrix basisColumnMatrix(basisDecompositionVector);
        // Can now compute the weights vector alpha
        Matrix Alpha = H.inverseLU() * basisColumnMatrix;
        // And finally update the stopping times - need to re-loop through all the paths
        for (size_t nSimul = 0; nSimul < NbSims; nSimul++) {
            // Computing the basis expansion with alpha coefficients
            double alphaBasisExpansion = 0.0;
            double payoffSimulation = std::max<double>(basketValues[nSimul]->GetValue(currentTimestep) - K, 0);
            for (size_t l=0; l < L; l++){
                alphaBasisExpansion += Alpha.data[l][0] * basisVectors[nSimul][l];
            }
            if (payoffSimulation >= alphaBasisExpansion){
                stoppingTimes[nSimul][i] = currentTimestep;
            }
            else {
                stoppingTimes[nSimul][i] = stoppingTimes[nSimul][i+1];
            }
        }
    }
    // Finally we can value the option
    vector<double> Payoffs(NbSims, 0.0);
    for (size_t nSimul = 0; nSimul < NbSims; nSimul++) {
        double stoppingTimeSimul = stoppingTimes[nSimul][1]; // at t=0, can't update the stopping times since all paths are the same => H matrix non invertible
        if (!UseControlVariate)
            Payoffs[nSimul] = exp(-Rate * stoppingTimeSimul) *
                              std::max<double>(basketValues[nSimul]->GetValue(stoppingTimeSimul) - K, 0);
        else {
            double ControlVariateLocalPayoff = 0.0;
            for (size_t d = 0; d < Dimension; d++) {
                ControlVariateLocalPayoff += Weights[d] * log(AllAssetPaths[nSimul][d]->GetValue(T));
            }
            Payoffs[nSimul] = TheoreticalPrice + exp(-Rate * stoppingTimeSimul) * (std::max<double>(
                    basketValues[nSimul]->GetValue(stoppingTimeSimul) - K, 0) - max<double>(
                    exp(ControlVariateLocalPayoff) - K, 0.0));
        }
    }
    double meanPrice = meanVector(Payoffs);
    double varPrice = varianceVector(Payoffs);
    end = clock();
    cout << "The price of the Bermudean Basket Call with MC is : " << meanPrice << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
    cout << "The variance of the Bermudean Basket Call with MC is : " << varPrice << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
    std::vector<double> results;
    results.reserve(3);
    results.push_back(meanPrice);
    results.push_back(varPrice);
    results.push_back(static_cast<double>(NbSims));
    return results;
};



BermudanBasketOption::~BermudanBasketOption()
{
}