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

double KakutaniSequence::Generate() {
    return 1.0;
}
