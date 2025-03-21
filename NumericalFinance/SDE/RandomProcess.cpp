#include "RandomProcess.h"

RandomProcess::RandomProcess(RandomGenerator* Gen, int dim)
{
	Generator = Gen;
	Dimension = dim;
}

SinglePath* RandomProcess::GetPath(int dimension)
{
	return Paths[dimension];
}
RandomProcess::~RandomProcess()
{
}
