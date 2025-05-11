#include "KakutaniSequence.h"

KakutaniSequence::KakutaniSequence(int nbSims, int dim, int length):
NbSims(nbSims), Dimension(dim), Length(length), localD(0), localN(0)
{
    countNbSim = 0;
    firstDPrimeNumbers = firstDPrimes(); /* Computed only once */
    pAdicObjects.resize(Dimension);
    for (int i = 0; i < Dimension; ++i)
        pAdicObjects[i] = new PAdic(firstDPrimeNumbers[i]);
    // Preparing the 3D vector
    Sequence.resize(nbSims);
    for (int sim = 0; sim < nbSims; ++sim) {
        Sequence[sim].resize(Length);
        for (int t = 0; t < Length; ++t) {
            Sequence[sim][t].resize(dim);
        }
    }
    createKakutaniSequence3D();
}

// Return the first d prime numbers
std::vector<int> KakutaniSequence::firstDPrimes() {
    std::vector<int> primes;
    int num = 2;
    while (primes.size() < Dimension) {
        bool is_prime = true;
        for (int i = 2; i <= sqrt(num); ++i)
            if (num % i == 0) {
                is_prime = false;
                break;
            }
        if (is_prime)
            primes.push_back(num);
        ++num;
    }
    return primes;
}

void KakutaniSequence::createKakutaniSequence3D() {
    /* Want to fill the 3d vector only once at instantiation to save time */
    std::vector<std::vector<std::vector<double>>> result(NbSims,
    std::vector<std::vector<double>>(Length, std::vector<double>(Dimension)));

    std::vector<double> x(Dimension), y(Dimension), xi(Dimension);

    // Compute y_i = 1/p + 1/p^2
    for (int i = 0; i < Dimension; ++i) {
        int p = firstDPrimeNumbers[i];
        x[i] = 1.0 / p;
        y[i] = 1.0 / p + 1.0 / (p * p);
        xi[i] = x[i];  // Starting value for each dimension
    }

    // Fill the 3D sequence
    for (int sim = 0; sim < NbSims; ++sim) {
        std::vector<double> xiSim = xi;  // Copy starting state for this simulation

        for (int t = 0; t < Length; ++t) {
            for (int d = 0; d < Dimension; ++d) {
                result[sim][t][d] = xiSim[d];
                xiSim[d] = pAdicObjects[d]->add(xiSim[d], y[d]);
            }
        }

        // For shifted index logic: advance by nbSteps to avoid overlaps
        for (int d = 0; d < Dimension; ++d) {
            for (int s = 0; s < Length; ++s)
                xi[d] = pAdicObjects[d]->add(xi[d], y[d]);
        }
    }
    Sequence = result;
}

double KakutaniSequence::Generate() {
    /* Once we're here we already have one a nbSIms *nbSteps*d matrix of RVs.
     * We just need to return them in the correct order */
    // Safety check: make sure we don't go out of bounds
    if (countNbSim >= Sequence.size()) {
        throw std::out_of_range("All simulations exhausted in KakutaniSequence::Generate().");
    }

    double output = Sequence[countNbSim][localN][localD];

    // Move to next dimension or timestep
    if (localD == Dimension - 1) {
        localD = 0;
        localN += 1;

        // If we've exhausted this simulation's time steps, move to next simulation
        if (localN == Length) {
            localN = 0;
            countNbSim += 1;
        }
    } else {
        localD += 1;
    }
    return output;
}
