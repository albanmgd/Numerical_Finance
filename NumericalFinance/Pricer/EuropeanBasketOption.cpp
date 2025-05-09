#include "EuropeanBasketOption.h"
#include "../Utils/basic_functions.h"

#include <memory>
#include "../Utils/Matrix.h"
#include "../SDE/BSEulerND.h"


EuropeanBasketOption::EuropeanBasketOption(
        size_t dim, double K, double T, double Rate, std::vector<double> Spots,
        std::vector<double> Vols, std::vector<double> Weights,
        std::vector<std::vector<double>> Correls, Normal* Gen) :
        BasketOption(dim, K, T, Rate, Spots, Vols,Weights,Correls, Gen)
{};

void EuropeanBasketOption::PriceCall(size_t NbSteps, size_t NbSims, bool UseAntithetic, bool UseControlVariate){
    cout << "Starting the MC Simulation ..." << endl;
    clock_t start, end;
    start = clock();

    BSEulerND TestScheme = BSEulerND(Generator, Dimension, Spots, Rate, Vols, &Correls);
    double TheoreticalPrice = TestScheme.PriceBasketCallOption(K, Weights, T, Correls);

    vector<double> Payoffs (NbSims, 0.0);

    for (size_t nSimul=0; nSimul < NbSims; nSimul++){
        double LocalPayoff = 0.0;
        double ControlVariateLocalPayoff = 0.0;
        TestScheme.Simulate(0, T, NbSteps, true);
        for (size_t d=0; d < Dimension; d++){
            LocalPayoff += Weights[d] * TestScheme.GetPath(d)->GetValue(T);
            ControlVariateLocalPayoff += Weights[d] * log(TestScheme.GetPath(d)->GetValue(T));
        }
        if (UseControlVariate == false)
            Payoffs[nSimul] = exp(-Rate * T) * max<double>(LocalPayoff - K, 0.0);
        else

            Payoffs[nSimul] =  exp(-Rate * T) * (max<double>(LocalPayoff - K, 0.0) - max<double>(exp(ControlVariateLocalPayoff) - K, 0.0)) + TheoreticalPrice;

    }
    end = clock();
    cout << "The price of the European Basket Call with MC is : " << meanVector(Payoffs) << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
    end = clock();
    cout << "The variance of the European Basket Call with MC is : " << varianceVector(Payoffs) << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
};

EuropeanBasketOption::~EuropeanBasketOption()
{
}
