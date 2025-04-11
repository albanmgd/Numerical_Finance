#include "BlackScholesND.h"

BlackScholesND::BlackScholesND(Normal* Gen, vector<double> spots, double rate, vector<double> vols) :
        RandomProcess(Gen, 1), Spots(spots), Rate(rate), Vols(vols)
{}

BlackScholesND::~BlackScholesND(){}

