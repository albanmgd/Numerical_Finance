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


void TestKakutaniSequence();
void TestVarianceReductionKakutaniSequence();
void TestLongstaffSchwarz();

int main()
{

//    TestKakutaniSequence();
//    TestVarianceReductionKakutaniSequence();
    TestLongstaffSchwarz();
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
    size_t nbSims = 1e4;
    std::vector<double> Spots = {100, 50, 60};
    std::vector<double> Vols = {0.10, 0.25, 0.16};
    double Rate = 0.05;
    std::vector<double> Weights = {0.10, 0.7, 0.2};
    std::vector<std::vector<double>> TestCorrelMatrix(dim, std::vector<double>(dim, 0.1));
    for (int i = 0; i < dim; ++i) {
        TestCorrelMatrix[i][i] = 1.0;
    }

    /* Params defining our pricing */
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);
    BSEulerND TestScheme = BSEulerND(NormBox, dim, Spots, Rate, Vols, &TestCorrelMatrix);
    int L = 3; /* Param controlling the expansion order for the basis functions */

    /* Begin by generating directly nbSim asset paths for the d assets */
    std::vector<std::vector<std::unique_ptr<SinglePath>>> AllAssetPaths(nbSims);
    for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
        AllAssetPaths[nSimul].resize(dim); /* resizing to hold d SinglePath */
        TestScheme.Simulate(0, T, nbSteps, false); /* Simulating the paths - enables us to reuse antithetic control variate*/
        for (size_t nbAsset = 0; nbAsset < dim; nbAsset++) {
            AllAssetPaths[nSimul][nbAsset] = std::make_unique<SinglePath>(*TestScheme.GetPath(nbAsset)); /* need to do copy of paths otherwise pointers get deleted every time we re-call the generate function*/
        }
    }
    // Backward induction for the stopping times
    std::vector<std::vector<double>> stoppingTimes(nbSims, std::vector<double>(dim, T)); /* initialization */

    // Keeping in memory the baskets paths to avoid recomputing them too
    std::vector<std::vector<double>> basketValues(nbSteps, std::vector<double>(nbSims, 0.0));
    // Now for each timestep we need to get the weights
    double delta_t = T / nbSteps;
    std::vector<std::vector<double>> basisVectors(nbSims, std::vector<double>(L, 0.0)); // used to construct phi

    for (int i=(nbSteps - 1); i >= 0; i--){
        double currentTimestep = i * delta_t; double previousTimestep = (i + 1) * delta_t; /* going backward */
        // Initializing the right hand part of the sum to optimize
        std::vector<double> basisDecompositionVector (L, 0.0);

        for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
            for (size_t j=0; j < dim; j++){
                basketValues[i][nSimul] += Weights[j] * AllAssetPaths[nSimul][j]->GetValue(currentTimestep);
            }
            double basketValue = basketValues[i][nSimul];
            double previousStoppingTime = stoppingTimes[nSimul][i];
            double mulFactor = exp(- Rate * (previousStoppingTime - currentTimestep)) *  std::max<double>(basketValue - K, 0);;
            /* Constructing the P vector */

            for (size_t l=0; l < L; l++){
                double basisScalar= pow(basketValue, l);
                basisVectors[nSimul][l] = basisScalar;
                basisDecompositionVector[l] += mulFactor * basisScalar;
            }
        }
        // Can now construct the Phi Matrix
        Matrix Phi = Matrix(basisVectors);
        Matrix H = Matrix(Phi.getTranspose() * Phi);
        // Some gymnastic with matrixes to build a column matrix
        std::vector<std::vector<double>> columnMatrix(L, std::vector<double>(1));
        for (size_t i = 0; i < L; ++i) {
            columnMatrix[i][0] = basisDecompositionVector[i];
        }
        Matrix basisColumnMatrix(columnMatrix);
        // Can now compute the weights vector alpha
        Matrix Alpha = H.inverseLU() * basisColumnMatrix;
        // And finally update the stopping times - need to reloop through all the paths
        for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
            // Computing the basis expansion with alpha coeffs
            double alphaBasisExpansion = 0.0; double payoffSimulation = std::max<double>(basketValues[i][nSimul] - K, 0);
            for (size_t l=0; l < L; l++){
                alphaBasisExpansion += Alpha.data[l][0] * basisVectors[nSimul][l];
            }
            if (payoffSimulation >= alphaBasisExpansion) stoppingTimes[nSimul][i] = currentTimestep;
            else stoppingTimes[nSimul][i] = stoppingTimes[nSimul][i+1];
        }
    }
    // Finally we can value the option
    double payoff = 0.0;
    for (size_t nSimul = 0; nSimul < nbSims; nSimul++) {
        double stoppingTimeSimul = stoppingTimes[nSimul][0];
        int indexStoppingTimeSimul = distance(stoppingTimes[nSimul].begin(),find(stoppingTimes[nSimul].begin(), stoppingTimes[nSimul].end(), stoppingTimeSimul));

        payoff += exp(- Rate * stoppingTimeSimul) * std::max<double>(basketValues[indexStoppingTimeSimul][nSimul] - K, 0);
    }
    cout << "The price of the Bermudean Basket Call is : " << payoff / nbSims << endl;
};