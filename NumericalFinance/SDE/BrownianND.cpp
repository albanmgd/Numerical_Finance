#include "BrownianND.h"
#include <vector>
#include <cmath>
#include <stdexcept>

BrownianND::BrownianND (RandomGenerator* Gen , int dim, vector <vector<double>>* Corr):
    RandomProcess(Gen, dim), Dimension(dim), CorrelationMatrix(Corr)
{
    CholeskyDecomposition = getCholeskyDecomposition(CorrelationMatrix); /* Computing it once only */
};

std::vector<std::vector<double>> BrownianND::getCholeskyDecomposition(std::vector<std::vector<double>>* Corr) {
    int n = Corr->size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = (*Corr)[i][j];
            for (int k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];

            if (i == j) {
                if (sum <= 0.0)
                    throw std::runtime_error("Matrix is not positive definite.");
                L[i][j] = std::sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return L;
}

void BrownianND::Simulate (double startTime, double endTime, size_t nbSteps){
    // Remove previous path
    if(Paths.size() > 0) delete Paths[0];
    Paths.clear();

/*    *//* Computing the Cholesky decomposition only once *//*
    std::vector<std::vector<double>> CholeskyMatrix = getCholeskyDecomposition(CorrelationMatrix);*/

    double dt = (endTime - startTime) / nbSteps;
    vector<vector<double>> MatrixCorrelatedBMs(Dimension, std::vector<double>(nbSteps, 0.0));

    /* We simulate the correlated brownian motions at every time step */
    for (size_t i = 0; i < nbSteps; ++i)
    {
        /* First need to generate the n random variables for a given timestamp */
        std::vector<double> NormalVariables(Dimension, 0.0);
        for (size_t k=0; k < Dimension; k++){
            NormalVariables[k] = Generator->Generate() * sqrt(dt);
        }

        /* Now we can compute the correlated Brownian motions - equivalent to matrix multiplication */
        for (size_t k=0; k < Dimension; k++){
            double CorrelatedBrownian = 0.0;
            /* Corresponds to matrix multiplication */
            for (size_t j=0; j < Dimension; j++){
                CorrelatedBrownian += CholeskyDecomposition[k][j] * NormalVariables[j];
            }
//            Path->AddValue(CorrelatedBrownian);
            MatrixCorrelatedBMs[k][i] = CorrelatedBrownian;
        }
    }
    // Want to store Paths as one Path = time series of Brownian motions for one asset
    for (size_t d=0; d < Dimension; d++){
        SinglePath* BrownianPath = new SinglePath(startTime, endTime, nbSteps);
        for (size_t i=0; i < nbSteps; i++){
            BrownianPath->AddValue(MatrixCorrelatedBMs[d][i]);
        }
        Paths.push_back(BrownianPath);
    }
}