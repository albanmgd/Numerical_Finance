#include "EuropeanBasketOption.h"
#include "../Utils/basic_functions.h"
#include "../Utils/Matrix.h"
#include "../SDE/BSEulerND.h"
#include <memory>

EuropeanBasketOption::EuropeanBasketOption(
        size_t Dim, double K, double T, double Rate, std::vector<double> Spots,
        std::vector<double> Vols, std::vector<double> Weights,
        std::vector<std::vector<double>> Correls, Normal* Gen) :
        BasketOption(Dim, K, T, Rate, Spots, Vols,Weights,Correls, Gen)
{};

std::vector<double> EuropeanBasketOption::PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic,
                                                    bool UseControlVariate){
    cout << "Starting the MC Simulation for the European Call..." << endl;
    clock_t start, end;
    start = clock();

    BSEulerND TestScheme = BSEulerND(Generator, Dimension, Spots, Rate, Vols, &Correls);
    double TheoreticalPrice = TestScheme.PriceBasketCallOption(K, Weights, T, Correls);

    vector<double> Payoffs (NbSims, 0.0);

    for (size_t nSimul=0; nSimul < NbSims; nSimul++){
        cout << "Done for sim: " << nSimul << endl;

        double LocalPayoff = 0.0;
        double ControlVariateLocalPayoff = 0.0;
        TestScheme.Simulate(0, T, NbSteps, UseAntithetic);
        for (size_t d=0; d < Dimension; d++){
            LocalPayoff += Weights[d] * TestScheme.GetPath(d)->GetValue(T);
            ControlVariateLocalPayoff += Weights[d] * log(TestScheme.GetPath(d)->GetValue(T));
        }
        if (!UseControlVariate)
            Payoffs[nSimul] = exp(-Rate * T) * max<double>(LocalPayoff - K, 0.0);
        else
            Payoffs[nSimul] =  exp(-Rate * T) * (max<double>(LocalPayoff - K, 0.0) - max<double>(exp(ControlVariateLocalPayoff) - K, 0.0)) + TheoreticalPrice;

    }
    end = clock();
    // storing results for output
    double meanPrice = meanVector(Payoffs);
    double variance = varianceVector(Payoffs);

    cout << "The price of the European Basket Call with MC is : " << meanPrice << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
    end = clock();
    cout << "The variance of the European Basket Call with MC is : " << variance << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;

    std::vector<double> results;
    results.reserve(3);
    results.push_back(meanPrice);
    results.push_back(variance);
    results.push_back(static_cast<double>(NbSims));
    return results;
};

EuropeanBasketOption::~EuropeanBasketOption()
{
}
