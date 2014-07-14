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
,surface_brightness(Data::get_instance().get_x_rays())
,model_image(Data::get_instance().get_ni(),
		vector<long double>(Data::get_instance().get_nj()))
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

	if(randomU() <= 0.5)
	{
		logH += source.perturb();

		if(source.get_size_of_diff() < source.get_num_components())
			update_surface_brightness();
		else
			calculate_surface_brightness();
		calculate_model_image();
	}
	else if(randomU() <= 0.5)
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
	else
	{
		logH += lens.perturb();

		if(lens.get_blobs_flag() &&
			(lens.get_size_of_diff() < lens.get_num_components()))
			update_rays();
		else
			shoot_rays();

		calculate_surface_brightness();
		calculate_model_image();
	}

	return logH;
}

double MyModel::logLikelihood() const
{
	double logL = 0.;
	const vector< vector<double> >& image =
				Data::get_instance().get_image();

	const vector< vector<double> >& sigma =
				Data::get_instance().get_sigma();

	double var;
	for(size_t i=0; i<image.size(); i++)
	{
		for(size_t j=0; j<image[i].size(); j++)
		{
			var = sigma[i][j]*sigma[i][j] + sigma0*sigma0
					+ sigma1*model_image[i][j];
			logL += -0.5*log(2.*M_PI*var) -
			0.5*pow(image[i][j] - model_image[i][j], 2)/var;
		}
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	out<<' '<<sigma0<<' '<<sigma1<<' ';
	lens.print(out);
	source.print(out);

	// Make an image of the source (uses the ray resolution)
	const vector< vector<long double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<long double> >& y = Data::get_instance().get_y_rays();
	for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
			out<<source.evaluate(x[i][j], y[i][j])<<' ';

	for(size_t i=0; i<model_image.size(); i++)
		for(size_t j=0; j<model_image[i].size(); j++)
			out<<model_image[i][j]<<' ';
}

string MyModel::description() const
{
	return string("lens parameters (b, q, rc, xc, yc, theta, shear, theta_shear), sigma0, sigma1, mock image");
}

void MyModel::shoot_rays()
{
	const vector< vector<long double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<long double> >& y = Data::get_instance().get_y_rays();

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

	// Blur using the PSF
	const PSF& psf = Data::get_instance().get_psf();

	if(Data::get_instance().use_fft())
		psf.blur_image2(surface_brightness);
	else
		psf.blur_image(surface_brightness);
}

void MyModel::update_surface_brightness()
{
	vector< vector<long double> >
		delta_surface_brightness(surface_brightness.size(),
			vector<long double>(surface_brightness[0].size(), 0.));

	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			delta_surface_brightness[i][j] = source.evaluate_diff
							(xs[i][j], ys[i][j]);
		}
	}

	// Blur using the PSF
	const PSF& psf = Data::get_instance().get_psf();

	if(Data::get_instance().use_fft())
		psf.blur_image2(delta_surface_brightness);
	else
		psf.blur_image(delta_surface_brightness);

	// Add the delta to the surface brightness
	for(size_t i=0; i<surface_brightness.size(); i++)
		for(size_t j=0; j<surface_brightness[i].size(); j++)
			surface_brightness[i][j] += delta_surface_brightness[i][j];
}


void MyModel::update_rays()
{
	const vector< vector<long double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<long double> >& y = Data::get_instance().get_y_rays();

	double ax, ay;
	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			lens.alpha_diff(x[i][j], y[i][j], ax, ay);
			xs[i][j] -= ax;
			ys[i][j] -= ay;
		}
	}
}

void MyModel::calculate_model_image()
{
	int resolution = Data::get_instance().get_resolution();

	model_image.assign(Data::get_instance().get_ni(),
		vector<long double>(Data::get_instance().get_nj(), 0.));

	int ii, jj;
	double coeff = pow(static_cast<double>(resolution), -2);
	for(size_t i=0; i<xs.size(); i++)
	{
		ii = i/resolution;
		for(size_t j=0; j<xs[i].size(); j++)
		{
			jj = j/resolution;
			model_image[ii][jj] += coeff*surface_brightness[i][j];
		}
	}
}


#include <fstream>
void MyModel::test()
{
	RandomNumberGenerator::initialise_instance();

	MyModel m;

	fstream fout("output.txt", ios::out);
	for(int i=0; i<1000; i++)
	{
		RandomNumberGenerator::get_instance().set_seed(i);
		m.fromPrior();
		m.perturb();
		m.print(fout); fout<<' ';
		m.shoot_rays();
		m.calculate_surface_brightness();
		m.calculate_model_image();
		m.print(fout); fout<<endl;
		cout<<(i+1)<<endl;
	}

	fout.close();
}

