#include "R2R1Function.h"
#include <exception>

R2R1Function::R2R1Function()
{}

R2R1Function::~R2R1Function()
{}

NullFunction::NullFunction() : R2R1Function()
{}

double NullFunction::operator()(double x, double t)
{
	return 0.0;
}

// Diffusion functions
BSActualization::BSActualization(double rate) : R2R1Function(), Rate(rate)
{}

double BSActualization::operator()(double x, double t)
{
	throw std::exception("Not implemented");
}

BSVariance::BSVariance(double sigma) : R2R1Function(), Sigma(sigma)
{}

double BSVariance::operator()(double x, double t)
{
	throw std::exception("Not implemented");
}

BSTrend::BSTrend(double rate) : R2R1Function(), Rate(rate)
{}

double BSTrend::operator()(double x, double t)
{
	throw std::exception("Not implemented");
}