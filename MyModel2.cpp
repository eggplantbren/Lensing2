#include "MyModel2.h"
#include "Data.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace Lensing2;
using namespace DNest3;

MyModel2::MyModel2()
:source(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
,lens(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
,xs(Data::get_instance().get_x_rays())
,ys(Data::get_instance().get_y_rays())
,surface_brightness(Data::get_instance().get_x_rays())
,model_image(Data::get_instance().get_ni(),
		vector<long double>(Data::get_instance().get_nj()))
{

}

void MyModel2::fromPrior()
{
	source.from_prior();
	lens.from_prior();

	sigma0 = tan(M_PI*(0.97*randomU() - 0.485));
	sigma1 = tan(M_PI*(0.97*randomU() - 0.485));
	sigma0 = exp(sigma0); sigma1 = exp(sigma1);

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();
}

double MyModel2::perturb()
{
	double logH = 0.;

	if(randomU() <= 0.5)
	{
		logH += source.perturb();

		calculate_surface_brightness();
		calculate_model_image();
	}
	else if(randomU() <= 0.5)
	{
		sigma0 = log(sigma0);
		sigma0 = (atan(sigma0)/M_PI + 0.485)/0.97;
		sigma0 += randh();
		wrap(sigma0, 0., 1.);
		sigma0 = tan(M_PI*(0.97*sigma0 - 0.485));
		sigma0 = exp(sigma0);

		sigma1 = log(sigma1);
		sigma1 = (atan(sigma1)/M_PI + 0.485)/0.97;
		sigma1 += randh();
		wrap(sigma1, 0., 1.);
		sigma1 = tan(M_PI*(0.97*sigma1 - 0.485));
		sigma1 = exp(sigma1);
	}
	else
	{
		logH += lens.perturb();

		if(lens.get_blobs_flag() &&
			(lens.get_size_of_diff() < lens.get_num_components())
			&& staleness2 < 10)
			update_rays();
		else
			shoot_rays();

		calculate_surface_brightness();
		calculate_model_image();
	}

	return logH;
}

double MyModel2::logLikelihood() const
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

void MyModel2::print(std::ostream& out) const
{
	out<<setprecision(5);
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

	out<<staleness1<<' '<<staleness2<<' ';
}

string MyModel2::description() const
{
	return string("lens parameters (b, q, rc, xc, yc, theta, shear, theta_shear), sigma0, sigma1, mock image");
}

void MyModel2::shoot_rays()
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

	staleness2 = 0;
}

void MyModel2::calculate_surface_brightness()
{
	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			surface_brightness[i][j] = source.evaluate(xs[i][j],
								ys[i][j]);
		}
	}

	if(Data::get_instance().psf_is_highres())
	{
		// Blur using the PSF
		const PSF& psf = Data::get_instance().get_psf();

		if(Data::get_instance().use_fft())
			psf.blur_image2(surface_brightness);
		else
			psf.blur_image(surface_brightness);
	}

	staleness1 = 0;
}

void MyModel2::update_rays()
{
	staleness2++;

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

void MyModel2::calculate_model_image()
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

	if(!Data::get_instance().psf_is_highres())
	{
		// Blur using the PSF
		const PSF& psf = Data::get_instance().get_psf();

		if(Data::get_instance().use_fft())
			psf.blur_image2(model_image);
		else
			psf.blur_image(model_image);
	}
}


#include <fstream>
void MyModel2::test()
{
	RandomNumberGenerator::initialise_instance();

	MyModel2 m;

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

