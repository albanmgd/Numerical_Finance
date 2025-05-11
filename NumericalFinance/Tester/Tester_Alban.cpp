#include <iostream>
#include <algorithm>
#include <ctime>

#include "../Utils/Matrix.h"
#include "../Utils/basic_functions.h"
#include "../Utils/CSVWriter.h"
#include "../RandomGenerator/Normal.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"
#include "../RandomGenerator/KakutaniSequence.h"
#include "../RandomGenerator/Normal.h"
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
    int testNbSims = 1e4;
    int testDim = 3; /* testing with d assets */
    int testN = 365;
    bool isAverageOk, isVarianceOk;
    KakutaniSequence* TestKakutaniSq = new KakutaniSequence(testNbSims, testDim, testN);
    UniformGenerator* UnifEuropean = new EcuyerCombined();
    // NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., TestKakutaniSq);
    NormalCLT* NormCLT = new NormalCLT(0., 1., TestKakutaniSq);

/*
    isAverageOk = NormCLT -> TestMean(testNbSims, 0.01);
    isVarianceOk = NormCLT -> TestVariance(testNbSims, 0.01);
    cout << "Mean OK: " << isAverageOk << " Variance OK:  " << isVarianceOk << std::endl;
*/

    // Exporting the generated nbs as a csv to check distrib with Python
    std::vector<std::vector<double>> values(testNbSims);

    for (size_t i = 0; i < testNbSims; ++i) {
        values[i] = { 0, 0, NormCLT->Generate() };  // Initialize each inner vector with one value
    }
    // Getting the csv
    std::string filename = "C:\\Users\\mager\\Downloads\\test_normal_clt.csv";
    WriteCSV(values, filename);
}

void TestClassImplementation() {
    int dim = 3;
    double T = 1.; // Maturity
    double K = 60;
    size_t nbSteps = 365;
    size_t nbSims = 1e4;
    vector<double> Spots = {100, 50, 60};
    vector<double> Vols = {0.10, 0.25, 0.16};
    double Rate = 0.05;
    vector<double> Weights = {0.10, 0.7, 0.2};

    vector<vector<double>> TestCorrelMatrix(dim, vector<double>(dim, 0.1));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    //UniformGenerator* UnifEuropean = new EcuyerCombined();
    UniformGenerator *UnifEuropean = new KakutaniSequence(nbSims, dim, nbSteps);

    NormalBoxMuller *NormEuropean = new NormalBoxMuller(0., 1., UnifEuropean);
    //NormalCLT* NormEuropean = new NormalCLT(0., 1., UnifEuropean);

    // Testing for the european option
    EuropeanBasketOption testEuropeanBasketOption(dim, K, T, Rate, Spots, Vols, Weights,
                                                  TestCorrelMatrix, NormEuropean);
    testEuropeanBasketOption.PriceCall(nbSteps, nbSims, false, false);

    // Testing for the bermudean option
    size_t L = 3;
    //UniformGenerator* UnifBermudean = new KakutaniSequence(nbSims, dim, nbSteps);
    UniformGenerator* UnifBermudean = new EcuyerCombined();

    NormalBoxMuller* NormBoxBermudean = new NormalBoxMuller(0., 1., UnifBermudean);
    BermudeanBasketOption testBermudeanBasketOption(dim, K, T, Rate, Spots, Vols, Weights,
                                                   TestCorrelMatrix, NormBoxBermudean, L);
    testBermudeanBasketOption.PriceCall(nbSteps, nbSims, false, false);
};