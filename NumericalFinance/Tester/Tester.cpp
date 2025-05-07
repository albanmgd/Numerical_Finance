#include <iostream>
#include <algorithm>
#include <ctime>

#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/FiniteSet.h"
#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/Poisson.h"
#include "../RandomGenerator/PAdic.h"
#include "../PDE/PDEGrid2D.h"
#include "../SDE/BlackScholes1D.h"
#include "../SDE/BSEuler1D.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"

void TestPDE();
void TestRandom();
void TestSDE();
void TestBrownianND();
void TestPAdic();

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
    cout << "p-adic decomposition of " << a << " and " << b << " yields: " << pAdicDecomposition.add(a, b);
}
