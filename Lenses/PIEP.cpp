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

void PIEP::alpha(double x, double y, double& ax, double& ay) const
{
	// Rotate and center
	double xx =  (x - xc)*cos_theta + (y - yc)*sin_theta;
	double yy = -(x - xc)*sin_theta + (y - yc)*cos_theta;

	double rsq = q*xx*xx + yy*yy/q + rc*rc;
	double r = sqrt(rsq);
	double alphax = b*q*xx/r;
	double alphay = b*yy/q/r;

	// Rotate back
	ax = alphax*cos_theta - alphay*sin_theta;
	ay = alphax*sin_theta + alphay*cos_theta;

	// Go into shear coordinate system
	xx =  x*cos_theta_shear + y*sin_theta_shear;
	yy = -x*sin_theta_shear + y*cos_theta_shear;

	// Calculate external shear
	alphax =  shear*xx;
	alphay = -shear*yy;

	// Add external shear
	ax += alphax*cos_theta_shear - alphay*sin_theta_shear;
	ay += alphax*sin_theta_shear + alphay*cos_theta_shear;
}

void PIEP::from_prior()
{
	b = exp(log(1E-3) + log(1E3)*randomU())*scale;
	q = exp(log(0.1) + log(10.)*randomU());
	rc = exp(log(1E-3) + log(1E3)*randomU())*scale;

	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(randomU() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(randomU() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	theta = M_PI*randomU();
	cos_theta = cos(theta); sin_theta = sin(theta);

	shear = 0.05*tan(M_PI*(randomU() - 0.5));
	theta_shear = 2.*M_PI*randomU();
	cos_theta_shear = cos(theta_shear); sin_theta_shear = sin(theta_shear);
}

double PIEP::perturb()
{
	double logH = 0.;

	int which = randInt(7);

	if(which == 0)
	{
		b = log(b/scale);
		b += log(1E3)*randh();
		b = mod(b - log(1E-3), log(1E3)) + log(1E-3);
		b = scale*exp(b);
	}
	else if(which == 1)
	{
		q = log(q);
		q += log(10.)*randh();
		q = mod(q - log(0.1), log(10.)) + log(0.1);
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
		logH -= -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH -= -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));

		xc += (x_max - x_min)*randh();
		yc += (y_max - y_min)*randh();

		xc = mod(xc - x_min, x_max - x_min) + x_min;
		yc = mod(yc - y_min, y_max - y_min) + y_min;

		logH += -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH += -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));
	}
	else if(which == 4)
	{
		theta += M_PI*randh();
		theta = mod(theta, 2.*M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}
	else if(which == 5)
	{
		shear = 0.5 + atan(shear/0.05)/M_PI;
		shear += randh();
		shear = mod(shear, 1.);
		shear = 0.05*tan(M_PI*(shear - 0.5));
	}
	else
	{
		theta_shear += 2.*M_PI*randh();
		theta_shear = mod(theta_shear, 2.*M_PI);
		cos_theta_shear = cos(theta_shear); sin_theta_shear = sin(theta_shear);
	}

	return logH;
}

void PIEP::print(ostream& out) const
{
	out<<b<<' '<<q<<' '<<rc<<' '<<xc<<' '<<yc<<' '<<theta<<' ';
	out<<shear<<' '<<theta_shear;
}

