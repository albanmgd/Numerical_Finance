#pragma once
#include "QuasiRandomGenerator.h"
#include "UniformGenerator.h"
#include "PAdic.h"

class KakutaniSequence : public UniformGenerator
{
protected:
    int Dimension; /* nb of assets */
    int Length; /* nb of timesteps */
    std::vector<int> firstDPrimeNumbers; /* computing them only once - used to generate the Kakutani sq */
    std::vector<std::vector<double>> Sequence;
    int localD; /* will be used with localN to have the function Generate to return only a double */
    int localN;

public:
    KakutaniSequence(int dim, int length);
    std::vector<int> firstDPrimes();
    std::vector<std::vector<double>> createKakutaniSequence();
    double Generate();
};