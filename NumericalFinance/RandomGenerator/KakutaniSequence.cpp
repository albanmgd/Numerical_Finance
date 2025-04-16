#pragma once
#include "KakutaniSequence.h"

KakutaniSequence::KakutaniSequence(PAdic* adicDecomp, int dim, int length):
Dimension(dim), Length(length), pAdicDecomp(adicDecomp)
{
    firstDPrimeNumbers = first_dprimes(); /* Computed only once */
}

// Return the first d prime numbers
std::vector<int> KakutaniSequence::first_dprimes() {
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

std::vector<std::vector<double>> KakutaniSequence::Generate() {
    std::vector<std::vector<double>> seq(Length, std::vector<double>(Dimension));
    std::vector<double> x(Dimension), y(Dimension);

    // Set x_i = 1/p_i, y_i = 1/p_i + 1/p_i^2
    for (int i = 0; i < Dimension; ++i) {
        int p = firstDPrimeNumbers[i];
        x[i] = 1.0 / p;
        y[i] = 1.0 / p + 1.0 / (p * p);
    }

    for (int t = 0; t < Length; ++t) {
        for (int i = 0; i < Dimension; ++i) {
            int p = firstDPrimeNumbers[i];
            double xi = x[i];
            for (int k = 0; k < t; ++k)
                xi = pAdicDecomp -> add(&xi, &y[i]);  // Apply T^t
            seq[t][i] = xi;
        }
    }
    return seq;
}
