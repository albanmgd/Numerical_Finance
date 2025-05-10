#include "BasketOption.h"

BasketOption::BasketOption(size_t Dim, double K, double T, double Rate, std::vector<double> Spots,
std::vector<double> Vols, std::vector<double> Weights,
std::vector<std::vector<double>> Correls, Normal* Gen):
Dimension(Dim), K(K), T(T), Rate(Rate), Spots(Spots), Vols(Vols), Weights(Weights), Correls(Correls), Generator(Gen){
    if ((Dim != Spots.size()) | (Dim != Vols.size()) | (Dim != Weights.size()) | (Dim != Correls.size()))
        std::runtime_error("Error in the dimensions of the input.");
};

BasketOption::~BasketOption()
{
}