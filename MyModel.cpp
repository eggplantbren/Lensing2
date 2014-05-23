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

	sigma0 = exp(log(1E-3) + log(1E6)*randomU());
	sigma1 = exp(log(1E-3) + log(1E6)*randomU());

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which = randInt(3);

	if(which == 0)
	{
		logH += source.perturb();

		calculate_surface_brightness();
		calculate_model_image();
	}
	else if(which == 1)
	{
		logH += lens.perturb();

		shoot_rays();
		calculate_surface_brightness();
		calculate_model_image();
	}
	else
	{
		sigma0 = log(sigma0);
		sigma0 += log(1E6)*randh();
		sigma0 = mod(sigma0 - log(1E-3), log(1E6)) + log(1E-3);
		sigma0 = exp(sigma0);

		sigma1 = log(sigma1);
		sigma1 += log(1E6)*randh();
		sigma1 = mod(sigma1 - log(1E-3), log(1E6)) + log(1E-3);
		sigma1 = exp(sigma1);
	}

	return logH;
}

double MyModel::logLikelihood() const
{
	double logL = 0.;
	const vector< vector<double> >& image =
				Data::get_instance().get_image();

	double var;
	for(size_t i=0; i<image.size(); i++)
	{
		for(size_t j=0; j<image[i].size(); j++)
		{
			var = sigma0*sigma0 + sigma1*model_image[i][j];
			logL += -0.5*log(2.*M_PI*var) -
			0.5*pow(image[i][j] - model_image[i][j], 2)/var;
		}
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	lens.print(out); out<<' '<<sigma0<<' '<<sigma1<<' ';
	for(size_t i=0; i<model_image.size(); i++)
		for(size_t j=0; j<model_image[i].size(); j++)
			out<<model_image[i][j]<<' ';
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
	// Get the psf from the data instance
	const PSF& psf = Data::get_instance().get_psf();

	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			surface_brightness[i][j] = source.evaluate(xs[i][j],
								ys[i][j]);
		}
	}

	psf.blur_image(surface_brightness);
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

