#include "R1R1Function.h"
#include <cmath>
#include <algorithm>

R1R1Function::R1R1Function()
{}

R1R1Function::~R1R1Function()
{}

VanillaTerminalCondition::VanillaTerminalCondition(double strike) : 
	R1R1Function(), Strike(strike)
{}

CallTerminalCondition::CallTerminalCondition(double strike) : VanillaTerminalCondition(strike)
{}

double CallTerminalCondition::operator()(double x)
{
	throw std::exception("Not implemented");
}

PutTerminalCondition::PutTerminalCondition(double strike) : VanillaTerminalCondition(strike)
{}

double PutTerminalCondition::operator()(double x)
{
	throw std::exception("Not implemented");
}

// Top and Bottom conditions

CallTopBoundary::CallTopBoundary(double sMax, double strike) : R1R1Function(), SMax(sMax), Strike(strike)
{}

double CallTopBoundary::operator()(double t)
{
	throw std::exception("Not implemented");
}

PutTopBoundary::PutTopBoundary(double sMax, double strike) : R1R1Function(), SMax(sMax), Strike(strike)
{}

double PutTopBoundary::operator()(double t)
{
	throw std::exception("Not implemented");
}

CallBottomBoundary::CallBottomBoundary(double sMin, double strike) : R1R1Function(), SMin(sMin), Strike(strike)
{}

double CallBottomBoundary::operator()(double t)
{
	throw std::exception("Not implemented");
}

PutBottomBoundary::PutBottomBoundary(double sMin, double strike) : R1R1Function(), SMin(sMin), Strike(strike)
{}

double PutBottomBoundary::operator()(double t)
{
	throw std::exception("Not implemented");
}