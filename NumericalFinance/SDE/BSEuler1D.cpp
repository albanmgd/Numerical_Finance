#include "BSEuler1D.h"

BSEuler1D::BSEuler1D(Normal* Gen, double spot, double rate, double vol) :
	BlackScholes1D(Gen, spot, rate, vol)
{
}

BSEuler1D::~BSEuler1D()
{
}

void BSEuler1D::Simulate(double startTime, double endTime, size_t nbSteps)
{
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->AddValue(Spot);
	double dt = (endTime - startTime) / nbSteps;
	double lastInserted = Spot;
	for (size_t i = 0; i < nbSteps; ++i)
	{
		double nextValue = lastInserted
			+ lastInserted * (Rate * dt + Vol * Generator->Generate() * sqrt(dt));

		Path->AddValue(nextValue);
		lastInserted = nextValue;
	}
	// Remove previous path
	if(Paths.size() > 0) delete Paths[0];
	Paths.clear();
	// Store new path
	Paths.push_back(Path);
}