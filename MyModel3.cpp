#include "MyModel3.h"
#include "Data.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

const double MyModel3::x_min = -1.;
const double MyModel3::x_max =  1.;
const double MyModel3::y_min = -1.;
const double MyModel3::y_max =  1.;
const double MyModel3::L = sqrt((x_max - x_min)*(y_max - y_min));

MyModel3::MyModel3()
:lens(x_min, x_max, y_min, y_max)
{

}

void MyModel3::fromPrior()
{
	lens.from_prior();
}

double MyModel3::perturb()
{
	double logH = 0.;

	logH += lens.perturb();

	return logH;
}

double MyModel3::logLikelihood() const
{
	double logL = 0.;

	// Data: positions of images
	double x[3] = {0., -0.220, -0.060};
	double y[3] = {0., -0.307, -0.139};
	double f[3] = {2.577, 1.977, 0.441};

	double xs[3], ys[3];
	double ax, ay;
	for(int i=0; i<3; i++)
	{
		lens.alpha(x[i], y[i], ax, ay);
		xs[i] = x[i] - ax;
		ys[i] = y[i] - ay;
	}

	double xs_mean = 0.;
	double ys_mean = 0.;
	for(int i=0; i<3; i++)
	{
		xs_mean += xs[i];
		ys_mean += ys[i];
	}
	xs_mean /= 3;
	ys_mean /= 3;

	double xsq = 0.;
	double ysq = 0.;
	for(int i=0; i<3; i++)
	{
		xsq += pow(xs[i] - xs_mean, 2);
		ysq += pow(ys[i] - ys_mean, 2);
	}
	xsq /= 3;
	ysq /= 3;

	double scatter = sqrt(sqrt(xsq)*sqrt(ysq));
	logL = -scatter/1E-3;

	// Add term for fluxes
	double flux1 = magnification(x[0], y[0]);
	double flux;
	for(int i=1; i<3; i++)
	{
		flux = magnification(x[i], y[i]);
		flux /= flux1;
		logL += -0.5*pow((flux - f[i]/f[0])/0.1, 2);
	}

	// Make sure no images are formed elsewhere
	double dist;
	bool check;
	for(double yy=1.; yy >= -1; yy -= 0.01)
	{
		for(double xx = -1; xx <= 1.; xx += 0.01)
		{
			check = true;
			for(int k=0; k<3; k++)
			{
				dist = sqrt(pow(xx - x[k], 2) + pow(yy - y[k], 2));
				if(dist < 0.01)
					check = false;
			}

			if(check)
			{
				lens.alpha(xx, yy, ax, ay);
				dist = sqrt(pow(xx - ax - xs_mean, 2) + pow(yy - ay - ys_mean, 2));
				if(dist < 0.01 && magnification(xx, yy) > 0.1)
					return -1E300;
			}
		}
	}

	return logL;
}

void MyModel3::print(std::ostream& out) const
{
	out<<setprecision(5);
	lens.print(out);

	// Data: positions of images
	double x[3] = {0., -0.220, -0.060};
	double y[3] = {0., -0.307, -0.139};
	double mag = 0.;
	for(int i=0; i<3; i++)
		mag += magnification(x[i], y[i]);
	out<<mag<<' ';
}

string MyModel3::description() const
{
	return string("");
}

double MyModel3::magnification(double x, double y) const
{
	double h = 1E-6;

	double J11, J12, J21, J22;
	double xs1, xs2, ys1, ys2, ax, ay;

	lens.alpha(x - h, y, ax, ay);
	xs1 = x - ax;
	ys1 = y - ay;
	lens.alpha(x + h, y, ax, ay);
	xs2 = x - ax;
	ys2 = y - ay;
	J11 = (xs2 - xs1)/(2.*h);
	J21 = (ys2 - ys1)/(2.*h);

	lens.alpha(x, y - h, ax, ay);
	xs1 = x - ax;
	ys1 = y - ay;
	lens.alpha(x, y + h, ax, ay);
	xs2 = x - ax;
	ys2 = y - ay;
	J12 = (xs2 - xs1)/(2.*h);
	J22 = (ys2 - ys1)/(2.*h);

	return 1./abs(J11*J22 - J12*J21);
}

