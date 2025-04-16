#include "BSEuler1D.h"

BSEuler1D::BSEuler1D(Normal* Gen, double spot, double rate, double vol) :
	BlackScholes1D(Gen, spot, rate, vol)
{
}

BSEuler1D::~BSEuler1D()
{
}

void BSEuler1D::Simulate(double startTime, double endTime, size_t nbSteps, bool antitheticRV)
{
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->AddValue(Spot);
	double dt = (endTime - startTime) / nbSteps;
	double lastInserted = Spot;
	SinglePath* AntitheticPath = nullptr; 
	double lastInsertedAntithetic = 0.0; // not sure ?

	if (antitheticRV){
		// we create a new antithetic path

		AntitheticPath = new SinglePath(startTime, endTime, nbSteps);
		AntitheticPath->AddValue(Spot);
		lastInsertedAntithetic = Spot;
	}

	for (size_t i = 0; i < nbSteps; ++i)
	{
		double nextValue = lastInserted
			+ lastInserted * (Rate * dt + Vol * Generator->Generate() * sqrt(dt));

		Path->AddValue(nextValue);
		lastInserted = nextValue;
		
		if (antitheticRV){
			double nextValueAntithetic = lastInsertedAntithetic 
			+ lastInsertedAntithetic * (Rate * dt - Vol * Generator->Generate() * sqrt(dt));
			AntitheticPath->AddValue(nextValueAntithetic);
			lastInsertedAntithetic = nextValueAntithetic;
		}

	}
	// Remove previous paths
	for (SinglePath* path : Paths)
		delete path;
	Paths.clear();

	// Store new paths
	Paths.push_back(Path);
	if (antitheticRV && AntitheticPath)
		Paths.push_back(AntitheticPath);
}