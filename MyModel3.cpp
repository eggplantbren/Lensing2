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

