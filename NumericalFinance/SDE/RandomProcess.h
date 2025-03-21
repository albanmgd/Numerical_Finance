#pragma once

#include "SinglePath.h"
#include "../RandomGenerator/Normal.h"

class RandomProcess
{
protected:
	RandomGenerator* Generator;
	vector<SinglePath*> Paths;
	int Dimension;

public:
	RandomProcess(RandomGenerator* Gen, int dim);
	virtual void Simulate(double startTime, double endTime, size_t nbSteps) = 0;
	SinglePath* GetPath(int dimension);
	~RandomProcess();
};

