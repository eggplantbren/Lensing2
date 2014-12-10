#include "MyModel3.h"
#include "Data.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

MyModel3::MyModel3()
:lens(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
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
	double x[4] = {-0.56, -0.65, -0.56, 0.2};
	double y[4] = {-0.2, 0., 0.2, 0.};

	// Trace into source plane
	double xs[4], ys[4];
	double ax, ay;
	for(int i=0; i<4; i++)
	{
		lens.alpha(x[i], y[i], ax, ay);
		xs[i] = x[i] - ax;
		ys[i] = y[i] - ay;
	}

	// Measure spread
	double meanx = 0.;
	double meany = 0.;
	for(int i=0; i<4; i++)
	{
		meanx += xs[i];
		meany += ys[i];
	}
	meanx /= 4.;
	meany /= 4.;

	double mean_sq_dev_x = 0.;
	double mean_sq_dev_y = 0.;
	for(int i=0; i<4; i++)
	{
		mean_sq_dev_x += pow(xs[i] - meanx, 2);
		mean_sq_dev_y += pow(ys[i] - meany, 2);
	}
	mean_sq_dev_x /= 4.;
	mean_sq_dev_y /= 4.;

	double sd_x = sqrt(mean_sq_dev_x);
	double sd_y = sqrt(mean_sq_dev_y);
	logL = -sqrt(sd_x*sd_y);

	return logL;
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

