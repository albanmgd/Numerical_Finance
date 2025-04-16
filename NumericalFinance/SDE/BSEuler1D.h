#pragma once
#include "BlackScholes1D.h"
class BSEuler1D : public BlackScholes1D
{
public:
	BSEuler1D(Normal* Gen, double spot, double rate, double vol);
	void Simulate(double startTime, double endTime, size_t nbSteps, bool antitheticRV);
	~BSEuler1D();
	// std::vector<SinglePath*> Paths; // not sure ?
};

