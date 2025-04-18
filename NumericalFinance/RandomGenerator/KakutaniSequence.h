#pragma once
#include "UniformGenerator.h"
#include "PAdic.h"

class KakutaniSequence : public UniformGenerator
{
protected:
    int Dimension; /* nb of assets */
    int Length; /* nb of timesteps */
    int countNbSim; /* counting the nb of times I have to generate a Kakutani sequence */
    std::vector<int> firstDPrimeNumbers; /* computing them only once - used to generate the Kakutani sq */
    std::vector<std::vector<double>> Sequence; /* Our n * d Kakutani sequence for a given sequence */
    int localD; /* will be used with localN to have the function Generate to return only a double */
    int localN;

public:
    KakutaniSequence(int dim, int length);
    std::vector<int> firstDPrimes();
    void createKakutaniSequence(int shiftIndex);
    double Generate();
};