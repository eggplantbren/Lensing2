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
const double MyModel3::dist = 0.02; // Distance to integrate fluxes

MyModel3::MyModel3()
:lens(x_min, x_max, y_min, y_max)
{

}

void MyModel3::fromPrior()
{
	lens.from_prior();

	x_source = x_min + (x_max - x_min)*randomU();
	y_source = y_min + (y_max - y_min)*randomU();

	width_source = exp(log(1E-3*L) + log(1E3)*randomU());
}

double MyModel3::perturb()
{
	double logH = 0.;

	int which = randInt(3);

	if(which == 0)
	{
	 	logH += lens.perturb();
	}
	else if(which == 1)
	{
		x_source += (x_max - x_min)*randh();
		y_source += (y_max - y_min)*randh();
		wrap(x_source, x_min, x_max);
		wrap(y_source, y_min, y_max);
	}
	else
	{
		width_source = log(width_source);
		width_source += 1E-3*randh();
		wrap(width_source, log(1E-3*L), log(L));
		width_source = exp(width_source);
	}

	return logH;
}

double MyModel3::logLikelihood() const
{
	double logL = 0.;

	// Data: positions of images
	double x[4] = {-0.56, -0.65, -0.56, 0.2};
	double y[4] = {-0.2, 0., 0.2, 0.};
	double f[4] = {1., 1.5, 1.2, 0.3};

	vector<double> flux(4);
	for(int i=0; i<4; i++)
	{
		flux[i] = flux_near(x[i], y[i]);
		if(i > 0)
			flux[i] /= flux[0];
		logL += -0.5*pow((f[i] - flux[i])/0.05, 2);
	}

	if(flux[0] == 0.)
		return -1E300;

	return logL;
}

double MyModel3::flux_near(double x, double y) const
{
	// u = log(r), du = 1/r dr, dr = r du = e^u du
	// dx dy = r dr dphi = e^{2u} du dphi

	vector<double> log_r(101);
	vector<double> phi(100);

	for(size_t i=0; i<log_r.size(); i++)
		log_r[i] = -10. + i*10./(log_r.size() - 1);

	for(size_t i=0; i<phi.size(); i++)
		phi[i] = (2.*M_PI*i)/100;

	double flux = 0.;
	double wsq = pow(width_source, 2);

	double xx, yy, ax, ay, xs, ys;
	for(size_t i=0; i<log_r.size(); i++)
	{
		for(size_t j=0; j<phi.size(); j++)
		{
			xx = x + exp(log_r[i])*cos(phi[j]);
			yy = y + exp(log_r[i])*sin(phi[j]);
			lens.alpha(xx, yy, ax, ay);
			xs = xx - ax;
			ys = yy - ay;

			if(pow(xs - x_source, 2) + pow(ys - y_source, 2) < wsq)
				flux += exp(2.*log_r[i]);
		}
	}

	return flux;
}

void MyModel3::print(std::ostream& out) const
{
	out<<setprecision(5);
	lens.print(out);
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

