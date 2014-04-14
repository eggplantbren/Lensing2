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
,lens(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
,xs(Data::get_instance().get_x_rays())
,ys(Data::get_instance().get_y_rays())
{

}

void MyModel::fromPrior()
{
	source.from_prior();
	lens.from_prior();

	shoot_rays();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which = randInt(2);

	if(which == 0)
		logH += source.perturb();
	else
	{
		logH += lens.perturb();
		shoot_rays();
	}

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

void MyModel::shoot_rays()
{

}


