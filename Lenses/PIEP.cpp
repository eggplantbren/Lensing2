#include "PIEP.h"
#include "RandomNumberGenerator.h"

#include <cmath>
#include <cassert>

using namespace std;
using namespace DNest3;
using namespace Lensing2;

PIEP::PIEP(double x_min, double x_max, double y_min, double y_max)
:x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
,scale(sqrt((x_max - x_min)*(y_max - y_min)))
{
	assert(x_max > x_min && y_max > y_min);
}

void PIEP::alpha(double x, double y, double& ax, double& ay)
{

}

void PIEP::from_prior()
{
	b = exp(log(1E-3) + log(1E3)*randomU())*scale;
	q = exp(log(0.1) + log(1E2)*randomU());
	rc = exp(log(1E-3) + log(1E3)*randomU())*scale;

	xc = x_min + (x_max - x_min)*randomU();
	yc = y_min + (y_max - y_min)*randomU();

	theta = 2.*M_PI*randomU();
	cos_theta = cos(theta); sin_theta = sin(theta);
}

double PIEP::perturb()
{
	return 0.;
}

void PIEP::print(ostream& out) const
{

}

