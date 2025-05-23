#include "SinglePath.h"
#include <stdexcept>

SinglePath::SinglePath(double start, double end, size_t nbSteps) : 
	StartTime(start), 
	EndTime(end), 
	NbSteps(nbSteps)
{
	if (nbSteps == 0)
	{
		throw std::runtime_error("Nb Steps is zero");
	}
	timeStep = (EndTime - StartTime) / nbSteps;
}

void SinglePath::AddValue(double val)
{
	Values.push_back(val);
	if (Times.size() == 0)
	{
		Times.push_back(StartTime);
	}
	else
	{
		Times.push_back(Times.back() + timeStep); // last time + time step
	}
}

double SinglePath::GetValue(double time)
{
    double precision = 1e-10; // used to handle rounding errors at the end of the series
	double result = 0.;
	if (time <= StartTime)
	{
		result = Values[0];
	}
	else if (time >= Times.back() - precision) // time greater than EndTime or between the last value and EndTime up to precision
	{
		result = Values.back();
	}
	else
	{
		for (size_t i = 1; i < Times.size(); ++i)
		{
			if (Times[i] < time)
				continue;
			else if (Times[i] == time)
			{
				result = Values[i];
				break;
			}
			else
			{
				double upperTime = Times[i];
				double lowerTime = Times[i - 1];
				
				double upperValue = Values[i];
				double lowerValue = Values[i - 1];

				result = lowerValue * ((upperTime - time) / (upperTime - lowerTime))
					+
					upperValue * ((time - lowerTime) / (upperTime - lowerTime));
				break;
			}
		}
	}

	return result;
}

vector<double>& SinglePath::GetValues()
{
	return Values;
}


vector<double>& SinglePath::GetTimes()
{
	return Times;
}

SinglePath::~SinglePath()
{
}