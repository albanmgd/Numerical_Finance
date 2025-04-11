#include "Exponential.h"
#include <stdexcept>
Exponential::Exponential(double inputLambda, UniformGenerator* unif) : ContinuousGenerator(1./inputLambda, 1./(inputLambda* inputLambda), unif)
{
	if (inputLambda <= 0)
		throw std::runtime_error("Lambda must be strictly positive for Exponential distribution");
	Lambda = inputLambda;
}

Exponential::~Exponential()
{}

ExponentialInverseDistribution::ExponentialInverseDistribution(double inputLambda, UniformGenerator* unif) : Exponential(inputLambda, unif)
{}
ExponentialRejectionSampling::ExponentialRejectionSampling(double inputLambda, UniformGenerator* unif) : Exponential(inputLambda, unif)
{}

double ExponentialInverseDistribution::Generate()
{
	return -log(Uniform->Generate()) / Lambda;
}

double ExponentialRejectionSampling::Generate()
{
	throw std::runtime_error("The rejection sampling method is not implemented yet for Exponential distribution");
}