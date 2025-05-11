#include "BlackScholes1D.h"

BlackScholes1D::BlackScholes1D(Normal* Gen, double spot, double rate, double vol) :
	RandomProcess(Gen, 1), Spot(spot), Rate(rate), Vol(vol)
{

}

BlackScholes1D::~BlackScholes1D()
{
}
