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

