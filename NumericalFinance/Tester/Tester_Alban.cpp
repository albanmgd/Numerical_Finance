#include <iostream>
#include <algorithm>
#include <ctime>

#include "../RandomGenerator/Normal.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"
#include "../RandomGenerator/PAdic.h"
#include "../RandomGenerator/KakutaniSequence.h"

void TestPAdic();
void TestKakutaniSequence();
void TestVarianceReductionKakutaniSequence();

int main()
{
//    TestRandom();
//    TestPDE();
//    TestSDE();
//    TestBrownianND();
//    TestBSEulerND();
//    TestPAdic();
//    TestKakutaniSequence();
    TestVarianceReductionKakutaniSequence();
}

void TestKakutaniSequence(){
    int testDim = 3; /* testing with d assets over 100 timesteps */
    int testN = 100;
    KakutaniSequence TestKakutaniSq = KakutaniSequence(testDim, testN);
    for (size_t n=0; n < testN; n++){
        for (size_t d=0; d < testDim; d++){
            std::cout << "For time step: " << n << " and dim: " << d << " generated RV is: " << TestKakutaniSq.Generate() << std::endl;
        }
    }
};

void TestVarianceReductionKakutaniSequence(){
    cout << "Starting the MC Simulation ..." << endl;
    clock_t start, end;
    start = clock();

    int dim = 3;
    double T = 1.; // Maturity
    double K = 65;
    size_t nbSteps = 365;
    size_t nbSims = 1e4;
    std::vector<double> Spots = {100, 50, 60};
    std::vector<double> Vols = {0.10, 0.25, 0.16};
    double Rate = 0.05;
    std::vector<double> Weights = {0.10, 0.7, 0.2};

    std::vector<std::vector<double>> TestCorrelMatrix(dim, std::vector<double>(dim, 0.1));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    UniformGenerator* Unif = new KakutaniSequence(dim, nbSteps);
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);

    BSEulerND TestScheme = BSEulerND(NormBox, dim, Spots, Rate, Vols, &TestCorrelMatrix);

    double Payoffs = 0.0;
    for (size_t nSimul=0; nSimul < nbSims; nSimul++){
        double LocalPayoff = 0.0;
        TestScheme.Simulate(0, T, nbSteps);
        for (size_t d=0; d < dim; d++){
            LocalPayoff += Weights[d] * TestScheme.GetPath(d)->GetValue(T);
        }
        Payoffs += std::max<double>(LocalPayoff - K, 0.0);
    }
    double Price = exp(-Rate * T) * Payoffs / nbSims;
    end = clock();
    cout << "The price of the European Basket Call is : " << Price << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
};