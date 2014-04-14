#include "MyModel.h"
#include "Data.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

MyModel::MyModel()
:source(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max(),
	1E-3, 1E3)
{

}

void MyModel::fromPrior()
{
	source.from_prior();
}

double MyModel::perturb()
{
	double logH = 0.;

	logH += source.perturb();

	return logH;
}

double MyModel::logLikelihood() const
{
	return 0.;
}

void MyModel::print(std::ostream& out) const
{
	out<<' ';
}

string MyModel::description() const
{
	return string("");
}

