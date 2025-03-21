#pragma once
#include <vector>

using namespace std;
class SinglePath
{
protected:
	vector<double> Values;
	vector<double> Times;
	double StartTime;
	double EndTime;
	size_t NbSteps;
	double timeStep;

public:
	SinglePath(double start, double end, size_t nbSteps);
	void AddValue(double val);
	double GetValue(double time);
	vector<double>& GetValues();
	vector<double>& GetTimes();
	~SinglePath();
};

