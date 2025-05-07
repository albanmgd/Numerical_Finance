#include <iostream>
#include <algorithm>
#include <ctime>
#include <memory>

#include "../Utils/Matrix.h"
#include "../RandomGenerator/Normal.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"
#include "../RandomGenerator/KakutaniSequence.h"
#include "../RandomGenerator/EcuyerCombined.h"

void TestBSEulerND();
void TestKakutaniSequence();
void TestVarianceReductionKakutaniSequence();
void TestLongstaffSchwarz();

int main()
{
    TestBSEulerND();
//    TestKakutaniSequence();
//    TestVarianceReductionKakutaniSequence();
    TestLongstaffSchwarz();
}

void TestBSEulerND(){
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

    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);

    BSEulerND TestScheme = BSEulerND(NormBox, dim, Spots, Rate, Vols, &TestCorrelMatrix);

    double Payoffs = 0.0;
    for (size_t nSimul=0; nSimul < nbSims; nSimul++){
        double LocalPayoff = 0.0;
        TestScheme.Simulate(0, T, nbSteps, false);
        for (size_t d=0; d < dim; d++){
            LocalPayoff += Weights[d] * TestScheme.GetPath(d)->GetValue(T);
        }
        Payoffs += std::max<double>(LocalPayoff - K, 0.0);
    }
    double Price = exp(-Rate * T) * Payoffs / nbSims;
    end = clock();
    cout << "The price of the European Basket Call is : " << Price << " found in "
         << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << endl;
}

void TestKakutaniSequence(){
    int testNbSims = 3;
    int testDim = 3; /* testing with d assets over 100 timesteps */
    int testN = 5;
    KakutaniSequence TestKakutaniSq = KakutaniSequence(testDim, testN);
    for (size_t NbSim=0; NbSim < testNbSims; NbSim ++) {
        for (size_t n = 0; n < testN; n++) {
            for (size_t d = 0; d < testDim; d++) {
                std::cout << "Simulation: " << NbSim << ", time step: " << n << " and dim: " << d << " generated RV is: "
                          << TestKakutaniSq.Generate() << std::endl;
            }
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

void TestLongstaffSchwarz(){
    /* Initial params to define our option */
    int dim = 3;
    double T = 1.; // Maturity
    double K = 65;
    size_t nbSteps = 365;
    size_t nbSims = 1e3;
    std::vector<double> Spots = {100, 50, 60};
    std::vector<double> Vols = {0.10, 0.25, 0.16};
    double Rate = 0.05;
    std::vector<double> Weights = {0.10, 0.7, 0.2};
    std::vector<std::vector<double>> TestCorrelMatrix(dim, std::vector<double>(dim, 0.1));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    /* Params defining our random nb generators */
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);
    BSEulerND TestScheme = BSEulerND(NormBox, dim, Spots, Rate, Vols, &TestCorrelMatrix);
    int L = 3; /* Param controlling the expansion order for the basis functions */

    /* Begin by generating directly nbSim asset paths for the d assets & computing the basket values associated */
    double delta_t = T / nbSteps;
    vector<vector<std::unique_ptr<SinglePath>>> AllAssetPaths(nbSims);
    std::vector<SinglePath*> basketValues(nbSims, nullptr);
    for (size_t nSimul=0; nSimul < nbSims; nSimul++) {
        AllAssetPaths[nSimul].resize(dim); /* resizing to hold d SinglePath */
        TestScheme.Simulate(0, T, nbSteps, false); /* Simulating the paths - enables us to reuse antithetic control variate*/
        for (size_t nbAsset=0; nbAsset < dim; nbAsset++) {
            AllAssetPaths[nSimul][nbAsset] = std::make_unique<SinglePath>(*TestScheme.GetPath(nbAsset)); /* need to do copy of paths otherwise pointers get deleted every time we re-call the generate function*/
        }
        /* Can finally compute the basket values for that simulation */
        SinglePath* basketValueSimulation = new SinglePath(0.0, T, nbSteps);
        for (size_t i=0; i<nbSteps; i++){
            double currentTimestep = i * delta_t;
            double basketValueSim = 0.0;
            for (size_t d=0; d<dim; d++){
                basketValueSim += Weights[d] * AllAssetPaths[nSimul][d]->GetValue(currentTimestep);
            }
            basketValueSimulation->AddValue(basketValueSim);
        }
        basketValues[nSimul] = basketValueSimulation;
    }
    // Backward induction for the stopping times
    vector<vector<double>> stoppingTimes(nbSims, std::vector<double>(nbSteps, T));

    // Looping through time, starting at T
    for (size_t i=(nbSteps - 1); i >= 1; i--){
        std::vector<std::vector<double>> basisVectors(nbSims, std::vector<double>(L, 0.0));
        double currentTimestep = i * delta_t; double previousTimestep = (i + 1) * delta_t; /* going backward */
        // Initializing the right hand part of the sum to optimize
        std::vector<double> basisDecompositionVector (L, 0.0);
        // Going through all the simulations for one time-step
        for (size_t nSimul=0; nSimul < nbSims; nSimul++) {
            double basketValue = basketValues[nSimul]->GetValue(currentTimestep);
            double previousStoppingTime = stoppingTimes[nSimul][i + 1];
            double mulFactor = exp(- Rate * (previousStoppingTime - currentTimestep)) *  std::max<double>(basketValue - K, 0);
            /* Constructing the P vector */
            for (size_t l=0; l < L; l++){
                double basisScalar= pow(basketValue, l);
                basisVectors[nSimul][l] = basisScalar;
                basisDecompositionVector[l] += mulFactor * basisScalar;
            }
        }
        cout << "At i: " <<i << endl;
        // Can now construct the Phi Matrix
        Matrix Phi = Matrix(basisVectors);
        Matrix H = Matrix(Phi.getTranspose() * Phi);
        // Some gymnastic with matrices to build a column matrix - L usually not too big so should be OK
        std::vector<std::vector<double>> columnMatrix(L, std::vector<double>(1));
        for (size_t l = 0; l < L; l++) {
            columnMatrix[l][0] = basisDecompositionVector[l];
        }
        Matrix basisColumnMatrix(columnMatrix);
        // Can now compute the weights vector alpha
        Matrix Alpha = H.inverseLU() * basisColumnMatrix;
        // And finally update the stopping times - need to re-loop through all the paths
        for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
            // Computing the basis expansion with alpha coefficients
            double alphaBasisExpansion = 0.0; double payoffSimulation = std::max<double>(basketValues[nSimul]->GetValue(currentTimestep) - K, 0);
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
    double payoff = 0.0;
    for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
        double stoppingTimeSimul = stoppingTimes[nSimul][0];
        payoff += exp(- Rate * stoppingTimeSimul) * std::max<double>(basketValues[nSimul]->GetValue(stoppingTimeSimul) - K, 0);
    }
    cout << "The price of the Bermudean Basket Call is : " << payoff / nbSims << endl;
};