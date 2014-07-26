#include "MultiBandModel.h"
#include "MultiBandData.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

MultiBandModel::MultiBandModel()
:source(MultiBandData::get_instance().get_images().size(),
	Blobby(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max(),
	1E-3, 1E3))
,lens(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
{

}

void MultiBandModel::fromPrior()
{
	for(size_t i=0; i<source.size(); i++)
		source[i].from_prior();
	lens.from_prior();
}

double MultiBandModel::perturb()
{
	double logH = 0.;

	return logH;
}

double MultiBandModel::logLikelihood() const
{
	double logL = 0.;

	return logL;
}

void MultiBandModel::print(std::ostream& out) const
{
	out<<' ';
}

string MultiBandModel::description() const
{
	return string("");
}

