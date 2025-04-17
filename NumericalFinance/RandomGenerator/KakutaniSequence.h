#include "RandomGenerator.h"
#include "PAdic.h"

class KakutaniSequence :  public RandomGenerator
{
protected:
    int Dimension;
    int Length;
    PAdic* pAdicDecomp;
    std::vector<int> firstDPrimeNumbers;

public:
    KakutaniSequence(PAdic* adicDecomp, int dim, int length);
    std::vector<int> first_dprimes();
    double Generate();
};