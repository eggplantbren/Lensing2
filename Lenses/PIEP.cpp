#include "PIEP.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"

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
	double logH = 0.;

	int which = randInt(5);

	if(which == 0)
	{
		b = log(b);
		b += log(1E3)*randh();
		b = mod(b - log(1E-3), log(1E3)) + log(1E-3);
		b = exp(b);
	}
	else if(which == 1)
	{
		q = log(q);
		q += log(1E2)*randh();
		q = mod(q - log(0.1), log(1E2)) + log(0.1);
		q = exp(q);
	}
	else if(which == 2)
	{
		rc = log(rc/scale);
		rc += log(1E3)*randh();
		rc = mod(rc - log(1E-3), log(1E3)) + log(1E-3);
		rc = scale*exp(rc);
	}
	else if(which == 3)
	{
		xc += (x_max - x_min)*randh();
		yc += (y_max - y_min)*randh();

		xc = mod(xc - x_min, x_max - x_min) + x_min;
		yc = mod(yc - y_min, y_max - y_min) + y_min;
	}
	else
	{
		theta += 2.*M_PI*randh();
		theta = mod(theta, 2.*M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}

	return logH;
}

void PIEP::print(ostream& out) const
{
	out<<b<<' '<<q<<' '<<rc<<' '<<xc<<' '<<yc<<' '<<theta<<' ';
}

