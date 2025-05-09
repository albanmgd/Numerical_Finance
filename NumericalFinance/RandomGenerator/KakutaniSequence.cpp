#include "KakutaniSequence.h"

KakutaniSequence::KakutaniSequence(int dim, int length):
Dimension(dim), Length(length), localD(0), localN(0)
{
    countNbSim = 0;
    firstDPrimeNumbers = firstDPrimes(); /* Computed only once */
    createKakutaniSequence(countNbSim);
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

void KakutaniSequence::KakutaniSequence::createKakutaniSequence(int shiftIndex) {
    std::vector<std::vector<double>> seq(Length, std::vector<double>(Dimension));
    std::vector<double> x(Dimension), y(Dimension);

    // Set x_i = 1/p_i, y_i = 1/p_i + 1/p_i^2
    for (int i = 0; i < Dimension; ++i) {
        int p = firstDPrimeNumbers[i];
        x[i] = 1.0 / p;
        y[i] = 1.0 / p + 1.0 / (p * p);
    }
    std::vector<PAdic*> pAdicObjects(Dimension);
    // Creating the vector of p-adic objects once outside the loop
    for (int i = 0; i < Dimension; ++i)
        pAdicObjects[i] = new PAdic(firstDPrimeNumbers[i]);

    for (int t = 0; t < Length; ++t) {
        for (int i = 0; i < Dimension; ++i) {
            PAdic* pAdicDecomp = pAdicObjects[i];
            double xi = x[i];
            int tShifted = t + shiftIndex;  // Apply the shift
            for (int k = 0; k < tShifted; k++)
                xi = pAdicDecomp->add(xi, y[i]);  // Apply T^{tShifted}
            seq[t][i] = xi;
        }
    }
    // avoid memory leaks
    for (int i = 0; i < Dimension; ++i)
        delete pAdicObjects[i];
    Sequence = seq;
}

double KakutaniSequence::Generate() {
    /* Once we're here we already have one a n*d matrix of RVs.
     * We just need to send return one of them in the correct order */
    if ((localD == Dimension - 1) and (localN == Length - 1)){ /* we have n timesteps on d dimensions but index starts at 0 ...*/
        /* Means that I've already returned the full n * d matrix generated for a simulation */
        localD = 0; localN = 0;
        countNbSim += 1;
        createKakutaniSequence(countNbSim);
    }

    double output = Sequence[localN][localD];
    if (localD == Dimension - 1){
        localD = 0;
        localN += 1;
    }
    else{
        localD += 1;
    }
    return output;
}
