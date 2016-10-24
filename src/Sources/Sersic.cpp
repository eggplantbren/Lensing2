#include "Sersic.h"
#include "RandomNumberGenerator.h"

using namespace std;
using namespace Lensing2;
using namespace DNest3;

Sersic::Sersic(double x_min, double x_max,
					double y_min, double y_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,image_size(sqrt((x_max - x_min)*(y_max - y_min)))
{

}

double Sersic::evaluate(double x, double y) const
{
	double f = 0.;

	// In a coordinate system centered on the galaxy and aligned with
	// the galaxy
	double xx, yy;

	xx =  (x - xc)*cos_theta + (y - yc)*sin_theta;
	yy = -(x - xc)*sin_theta + (y - yc)*cos_theta;

	double rsq = q*xx*xx + yy*yy/q;
	f = Ie*exp(-bm*(pow(rsq*one_over_Re_squared, 0.5/m) - 1.));

	return f;
}


void Sersic::from_prior()
{
	// log(Ie) ~ Cauchy(0, 1)T(-21, 21)
	Ie = exp(tan(M_PI*(0.97*randomU() - 0.485)));

	// log(Re) ~ Uniform(log(1E-3*image_size), log(image_size))
	Re = exp(log(1E-3*image_size) + log(1E3)*randomU());
	one_over_Re_squared = 1./(Re*Re);

	m = 0.2 + 9.8*randomU();
	bm = 2.*m - 0.324;

	// xc and yc ~ Cauchy(center of image, 0.1*image_size)T(inside image)
	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(randomU() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(randomU() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	q = 0.1 + 0.9*randomU();
	theta = M_PI*randomU();
	cos_theta = cos(theta); sin_theta = sin(theta);
}

double Sersic::perturb()
{
	double logH = 0.;

	// randh gives you 10^(1.5 - 6*rand)*randn
	// wrap(x, a, b) replaces x with ((x-a) mod b) + a

	// Choose which parameter to modify
	int which = randInt(6); // Either 0, 1, 2, 3, 4, or 5.

	if(which == 0)
	{
		// Transforms from the crazy Ie prior to U(0, 1)
		// Then does proposal
		// Then transforms back
		Ie = log(Ie);
		Ie = (atan(Ie)/M_PI + 0.485)/0.97;
		Ie += randh();
		wrap(Ie, 0., 1.);
		Ie = tan(M_PI*(0.97*Ie - 0.485));
		Ie = exp(Ie);
	}
	if(which == 1)
	{
		m += 9.8*randh();
		wrap(m, 0.2, 10.);
		bm = 2.*m - 0.324;
	}
	if(which == 2)
	{
		// Transforms the Re value to log(Re) which has a
		// Uniform prior, then does proposal, then transforms back
		Re = log(Re);
		Re += log(1E3)*randh();
		wrap(Re, log(1E-3*image_size), log(image_size));
		Re = exp(Re);
		one_over_Re_squared = 1./(Re*Re);	
	}
	if(which == 3)
	{
		// Proposal respects the cauchy prior 
		logH -= -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH -= -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));

		xc += image_size*randh();
		yc += image_size*randh();
		wrap(xc, x_min, x_max);
		wrap(yc, y_min, y_max);

		logH += -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH += -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));
	}
	if(which == 4)
	{
 		q += 0.9*randh();
		wrap(q, 0.1, 1.);
	}
	
	if(which == 5)
	{
 		theta += M_PI*randh();
		wrap(theta, 0.0 , M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}


	return logH;
}

void Sersic::print(ostream& out) const
{
	out<<Ie<<' '<<m<<' '<<Re<<' '<<xc<<' '<<yc<<' '<<q<<' '<<theta<<' ';
}

