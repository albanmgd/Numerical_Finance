#include <iostream>
#include <algorithm>
#include <ctime>
#include <memory>

#include "../Utils/Matrix.h"
#include "../Utils/basic_functions.h"
#include "../RandomGenerator/Normal.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"
#include "../RandomGenerator/KakutaniSequence.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../Pricer/EuropeanBasketOption.h"
#include "../Pricer/BermudeanBasketOption.h"

void TestPAdic();
void TestMeanKakutani();
void TestClassImplementation();
void TestKakutaniSequence();
void TestVarianceReductionKakutaniSequence();

int main()
{
//   TestPAdic();
//   TestKakutaniSequence();
//    TestMeanKakutani();
    TestClassImplementation();
//    TestVarianceReductionKakutaniSequence();
//    TestLongstaffSchwarz();
}

void TestPAdic(){
    /* Testing the example of the course */
    double a = 0.123333333;
    double b = 0.412777777;
    PAdic pAdicDecomposition = PAdic(10);
    cout << "p-adic decomposition of " << a << " and " << b << " yields: " << pAdicDecomposition.add(a, b) << std::endl;;
}

void TestKakutaniSequence(){
    int testNbSims = 3;
    int testDim = 3; /* testing with d assets */
    int testN = 5;
    KakutaniSequence TestKakutaniSq = KakutaniSequence(testNbSims, testDim, testN);
    for (size_t NbSim=0; NbSim < testNbSims; NbSim ++) {
        for (size_t n = 0; n < testN; n++) {
            for (size_t d = 0; d < testDim; d++) {
                std::cout << "Simulation: " << NbSim << ", time step: " << n << " and dim: " << d << " generated RV is: "
                          << TestKakutaniSq.Generate() << std::endl;
            }
        }
    }
};

void TestMeanKakutani() {
    int testNbSims = 1e3;
    int testDim = 3; /* testing with d assets */
    int testN = 365;
    bool isAverageOk, isVarianceOk;
    KakutaniSequence TestKakutaniSq = KakutaniSequence(testNbSims, testDim, testN);
    isAverageOk = TestKakutaniSq.TestMean(testNbSims, 0.01);
    isVarianceOk = TestKakutaniSq.TestVariance(testNbSims, 0.01);
    cout << "Mean OK: " << isAverageOk << " Variance OK:  " << isVarianceOk << std::endl;;

}

void TestClassImplementation(){
    int dim = 3;
    double T = 1.; // Maturity
    double K = 60;
    size_t nbSteps = 365;
    size_t nbSims = 1e3;
    vector<double> Spots = {100, 50, 60};
    vector<double> Vols = {0.10, 0.25, 0.16};
    double Rate = 0.05;
    vector<double> Weights = {0.10, 0.7, 0.2};

    vector<vector<double>> TestCorrelMatrix(dim, vector<double>(dim, 0.1));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    //UniformGenerator* UnifEuropean = new EcuyerCombined();
    UniformGenerator* UnifEuropean = new KakutaniSequence(nbSims, dim, nbSteps);
    NormalBoxMuller* NormBoxEuropean = new NormalBoxMuller(0., 1., UnifEuropean);

    // Testing for the european option
    EuropeanBasketOption testEuropeanBasketOption(dim, K, T, Rate, Spots, Vols, Weights,
                                                  TestCorrelMatrix, NormBoxEuropean);
    testEuropeanBasketOption.PriceCall(nbSteps, nbSims, true, false);

    // Testing for the bermudean option
    size_t L = 3;
    UniformGenerator* Unif = new KakutaniSequence(nbSims, dim, nbSteps);
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);
    BermudeanBasketOption testBermudeanBasketOption(dim, K, T, Rate, Spots, Vols, Weights,
                                                   TestCorrelMatrix, NormBox, L);
    testBermudeanBasketOption.PriceCall(nbSteps, nbSims, false, false);
};