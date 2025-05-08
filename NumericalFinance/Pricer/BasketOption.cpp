#include "BasketOption.h"

BasketOption::BasketOption(size_t dim, double K, double T, double Rate, std::vector<std::vector<double>>* Spots,
std::vector<std::vector<double>>* Vols, std::vector<std::vector<double>>* Weights,
std::vector<std::vector<double>>* Correls):
Dimension(dim), K(K), T(T), Rate(Rate), Spots(Spots), Vols(Vols), Weights(Weights), Correls(Correls){

    if ((dim != Spots->size()) | (dim != Vols->size()) | (dim != Weights->size()) | (dim != Correls->size()))
        std::runtime_error("Error in the dimensions of the input.");
};

BasketOption::~BasketOption()
{
}