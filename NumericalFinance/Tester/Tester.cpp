#include <iostream>
#include <algorithm>
#include <ctime>

#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/FiniteSet.h"
#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/PAdic.h"
#include "../PDE/PDEGrid2D.h"
#include "../SDE/BlackScholes1D.h"
#include "../SDE/BSEuler1D.h"
#include "../SDE/BrownianND.h"
#include "../Pricer/EuropeanBasketOption.h"
#include "../Pricer/BermudeanBasketOption.h"
#include "../Utils/Matrix.h"
#include "../Utils/basic_functions.h"
#include "../Utils/CSVWriter.h"
#include "../RandomGenerator/KakutaniSequence.h"

void TestPDE();
void TestRandom();
void TestSDE();
void TestBrownianND();
void TestPAdic();
void TestMeanKakutani();
void TestClassImplementation();
void TestKakutaniSequence();
void TestVarianceReductionKakutaniSequence();
/*int main()
{
//    TestRandom();
//    TestPDE();
//    TestSDE();
//    TestBrownianND();
    TestPAdic();
}*/

void TestRandom()
{
    myLong Multiplier, Increment, Modulus, Seed;

    Multiplier = 40014;
    Increment = 0;
    Modulus = 2147483563;
    Seed = 1;

    //X = new LinearCongruential(Multiplier, Increment, Modulus, Seed);
    //X = new EcuyerCombined();
    UniformGenerator* Unif = new EcuyerCombined();
    std::vector<double> probs(3);
    probs[0] = 0.4;
    probs[1] = 0.5;
    probs[2] = 0.1;
    RandomGenerator* Finit = new FiniteSet(probs, Unif);
    /*
        for (size_t i = 0; i < 10; ++i)
        {
            std::cout << Finit->Generate() << std::endl;
        }

        ExponentialInverseDistribution* Exp = new ExponentialInverseDistribution(0.5, Unif);
        */
    bool isAverageOk, isVarianceOk;
    //isAverageOk = Exp->TestMean(100000, 0.01);
    //isVarianceOk = Exp->TestVariance(100000, 0.1);

    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);
    isAverageOk = NormBox->TestMean(100000, 0.01);
    isVarianceOk = NormBox->TestVariance(100000, 0.1);
}

void TestPDE()
{
    /*
        Evaluate a Call option in the Black-Scholes framework:
        Expected Price : 6.80495771
    */
    double Spot = 100;
    double Strike = 100;
    double Maturity = 1;
    double Rate = 0.05;
    double Volatility = 0.1;

    double SMin = 0;
    double SMax = 1000;

    R2R1Function* VarianceFunction = new BSVariance(Volatility);
    R2R1Function* TrendFunction = new BSTrend(Rate);
    R2R1Function* ActualizationFunction = new BSActualization(Rate);
    R2R1Function* SourceTermFunction = new NullFunction();

    R1R1Function* TopBoundaryFunction = new CallTopBoundary(SMax, Strike);
    R1R1Function* BottomBoundaryFunction = new CallBottomBoundary(SMin, Strike);
    R1R1Function* RightBoundaryFunction = new CallTerminalCondition(Strike);

    PDEGrid2DExplicit BlackScholesGrid(Maturity, SMin, SMax, 1000, 1.,
        VarianceFunction, TrendFunction, ActualizationFunction, SourceTermFunction,
        TopBoundaryFunction, BottomBoundaryFunction, RightBoundaryFunction);

    BlackScholesGrid.FillNodes();
    double priceAtZero = BlackScholesGrid.GetTimeZeroNodeValue(Spot);

    std::cout << "Price = " << priceAtZero << std::endl;
}

void TestSDE()
{
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);

    double spot = 100.;
    double rate = 0.05;
    double vol = 0.1;
    // Build the Euler scheme for 1D black scholes
    BSEuler1D Scheme = BSEuler1D(NormBox, spot, rate, vol);
    
    /* Compute the price of an option */

    double T = 1.; // Maturity
    double K = 100; // Strike
    size_t nbSteps = 365; // Number of steps from 0 to T for one path
    size_t nbSimul = 100000.; // Number of Monte Carlo simulations

    double Price = 0.; // Price of the Call option

    for (size_t i = 0; i < nbSimul; ++i)
    {
        // Simulate the current path
        Scheme.Simulate(0., T, nbSteps);
        // Get the last value of ST on the current path
        double ST = Scheme.GetPath(0)->GetValue(T);
        // Compute the payoff
        double payoff = std::max<double>(ST - K, 0);
        // Compute and store average
        Price += payoff / nbSimul;
    }
    // Discount the price
    Price *= exp(-rate * T);
    // Print the price
    cout << "The price of BS Call is : " << Price << endl;

    // Example of path
    std::cout << "------------- Example of path ------------------" << endl;
    // Simulate the process once, to test the path generation
    Scheme.Simulate(0., 1., 100);
    // Get the current path
    SinglePath* path = Scheme.GetPath(0);
    // Read times and values from the path
    vector<double>& pathTimes = path->GetTimes();
    vector<double>& pathValues = path->GetValues();
    // Print the path
    for (size_t i = 0; i < pathTimes.size(); ++i)
    {
        std::cout << pathTimes[i] << " " << pathValues[i] << std::endl;
    }

    string input = "";
    std::cin >> input;
}

void TestBrownianND(){
    int dim = 3;
    double T = 1.; // Maturity
    size_t nbSteps = 365;

    std::vector<std::vector<double>> TestCorrelMatrix(dim, std::vector<double>(dim, 0.0));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);

    BrownianND TestBrownianND = BrownianND(NormBox, dim, &TestCorrelMatrix);
    TestBrownianND.Simulate(0, T, nbSteps);
};

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
