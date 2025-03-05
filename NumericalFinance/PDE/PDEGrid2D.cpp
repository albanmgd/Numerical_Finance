#include "PDEGrid2D.h"
#include <exception>

PDEGrid2D::PDEGrid2D(
	double Maturity,
	double MinUnderlyingValue,
	double MaxUnderlyingValue,
	int NbTimeSteps,
	double StepForUnderlying,
	R2R1Function* VarianceFunction,
	R2R1Function* TrendFunction,
	R2R1Function* ActualizationFunction,
	R2R1Function* SourceTermFunction,
	R1R1Function* TopBoundaryFunction,
	R1R1Function* BottomBoundaryFunction,
	R1R1Function* RightBoundaryFunction
)
{
	throw std::exception("Not implemented");
}

PDEGrid2D::~PDEGrid2D()
{
}

void PDEGrid2D::FillRightBoundary()
{
	throw std::exception("Not implemented");
}

void PDEGrid2D::FillTopAndBottomBoundary()
{
	throw std::exception("Not implemented");
}

void PDEGrid2D::FillNodes()
{
	throw std::exception("Not implemented");
}

double PDEGrid2D::GetTimeZeroNodeValue(double spot)
{
	throw std::exception("Not implemented");
}

PDEGrid2DExplicit::PDEGrid2DExplicit(
	double Maturity,
	double MinUnderlyingValue,
	double MaxUnderlyingValue,
	int NbTimeSteps,
	double StepForUnderlying,
	R2R1Function* VarianceFunction,
	R2R1Function* TrendFunction,
	R2R1Function* ActualizationFunction,
	R2R1Function* SourceTermFunction,
	R1R1Function* TopBoundaryFunction,
	R1R1Function* BottomBoundaryFunction,
	R1R1Function* RightBoundaryFunction
) : PDEGrid2D(
	Maturity,
	MinUnderlyingValue,
	MaxUnderlyingValue,
	NbTimeSteps,
	StepForUnderlying,
	VarianceFunction,
	TrendFunction,
	ActualizationFunction,
	SourceTermFunction,
	TopBoundaryFunction,
	BottomBoundaryFunction,
	RightBoundaryFunction
)
{}

void PDEGrid2DExplicit::FillNodes()
{
	throw std::exception("Not implemented");
}