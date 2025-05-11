#pragma once
#include "UniformGenerator.h"
#include "PAdic.h"

class KakutaniSequence : public UniformGenerator
{
protected:
    int NbSims;
    int Dimension; /* nb of assets */
    int Length; /* nb of timesteps */
    int countNbSim; /* counting the nb of times I have to generate a Kakutani sequence */
    std::vector<int> firstDPrimeNumbers; /* computing them only once - used to generate the Kakutani sq */
    std::vector<std::vector<std::vector<double>>> Sequence; /* Our n * d Kakutani sequence for a given sequence */
    int localD; /* will be used with localN to have the function Generate to return only a double */
    int localN;
    std::vector<PAdic*> pAdicObjects; /* Avoid having to re-instantiate them every time we regenerate a sq*/

public:
    KakutaniSequence(int nbSims, int dim, int length);
    std::vector<int> firstDPrimes();
    void createKakutaniSequence3D();
    /*void updateKakutaniSequence();*/
    double Generate();
};