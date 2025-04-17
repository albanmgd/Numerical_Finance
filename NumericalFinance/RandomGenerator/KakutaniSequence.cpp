#include "KakutaniSequence.h"

KakutaniSequence::KakutaniSequence(int dim, int length):
Dimension(dim), Length(length), localD(0), localN(0)
{
    firstDPrimeNumbers = firstDPrimes(); /* Computed only once */
    Sequence = createKakutaniSequence();
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

std::vector<std::vector<double>> KakutaniSequence::createKakutaniSequence() {
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
            PAdic* pAdicDecomp = new PAdic(p);
            double xi = x[i];
            for (int k = 0; k < t; ++k)
                xi = pAdicDecomp -> add(&xi, &y[i]);  // Apply T^t
            seq[t][i] = xi;
        }
    }
    return seq;
}

double KakutaniSequence::Generate() {
    /* Once we're here we have our n sequences of d random numbers already computed. We just need to send return one of them in the
     * correct order */
    double output;
    if ((localD == Dimension) and (localN == Length)){
        throw std::runtime_error("Please increase the dimensions of the sequence. All generated numbers have been returned.");
    }
    else{
        output = Sequence[localD][localN];
        if (localD == Dimension){
            localD = 0;
            localN += 1;
        }
        else{
            localD += 1;
        }
    }
    return output;
}
