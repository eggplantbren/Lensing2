#include "MultiBandModel.h"
#include "MultiBandData.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

MultiBandModel::MultiBandModel()
:source(MultiBandData::get_instance().get_images().size())
{

}

void MultiBandModel::fromPrior()
{

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

