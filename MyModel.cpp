#include "MyModel.h"
#include "Data.h"
#include "Constants.h"

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
,surface_brightness(Data::get_instance().get_x_rays())
,model_image(Data::get_instance().get_ni(),
		vector<double>(Data::get_instance().get_nj()))
{

}

void MyModel::fromPrior()
{
	source.from_prior();
	lens.from_prior();

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which = randInt(2);

	if(which == 0)
		logH += source.perturb();
	else
		logH += lens.perturb();

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();

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
	const vector< vector<double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<double> >& y = Data::get_instance().get_y_rays();

	double ax, ay;
	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			lens.alpha(x[i][j], y[i][j], ax, ay);
			xs[i][j] = x[i][j] - ax;
			ys[i][j] = y[i][j] - ay;
		}
	}
}

void MyModel::calculate_surface_brightness()
{
	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			surface_brightness[i][j] = source.evaluate(xs[i][j],
								ys[i][j]);
		}
	}
}

void MyModel::calculate_model_image()
{
	model_image.assign(Data::get_instance().get_ni(),
		vector<double>(Data::get_instance().get_nj(), 0.));

	int ii, jj;
	double coeff = pow(static_cast<double>(Constants::resolution), -2);
	for(size_t i=0; i<xs.size(); i++)
	{
		ii = i/Constants::resolution;
		for(size_t j=0; j<xs[i].size(); j++)
		{
			jj = j/Constants::resolution;
			model_image[ii][jj] += coeff*surface_brightness[i][j];
		}
	}
}

