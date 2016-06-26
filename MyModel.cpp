#include "MyModel.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest4;

namespace Lensing2
{

MyModel::MyModel()
:source(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
,lens(Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
	Data::get_instance().get_y_min(), Data::get_instance().get_y_max())
,xs(Data::get_instance().get_x_rays())
,ys(Data::get_instance().get_y_rays())
,surface_brightness(Data::get_instance().get_x_rays())
,model_image(Data::get_instance().get_ni(),
		vector<double>(Data::get_instance().get_nj()))
{

}

void MyModel::from_prior(RNG& rng)
{
	source.from_prior(rng);
	lens.from_prior(rng);

	sigma0 = tan(M_PI*(0.97*rng.rand() - 0.485));
	sigma1 = tan(M_PI*(0.97*rng.rand() - 0.485));
	sigma0 = exp(sigma0); sigma1 = exp(sigma1);

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();
}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;

	if(rng.rand() <= 0.5)
	{
		logH += source.perturb(rng);

		calculate_surface_brightness(source.get_blobs().get_removed().size() == 0);
		calculate_model_image();
	}
	else if(rng.rand() <= 0.5)
	{
		sigma0 = log(sigma0);
		sigma0 = (atan(sigma0)/M_PI + 0.485)/0.97;
		sigma0 += rng.randh();
		wrap(sigma0, 0., 1.);
		sigma0 = tan(M_PI*(0.97*sigma0 - 0.485));
		sigma0 = exp(sigma0);

		sigma1 = log(sigma1);
		sigma1 = (atan(sigma1)/M_PI + 0.485)/0.97;
		sigma1 += rng.randh();
		wrap(sigma1, 0., 1.);
		sigma1 = tan(M_PI*(0.97*sigma1 - 0.485));
		sigma1 = exp(sigma1);
	}
	else
	{
		logH += lens.perturb(rng);

		shoot_rays(lens.get_blobs_flag() &&
                   lens.get_blobs().get_removed().size() == 0);
		calculate_surface_brightness();
		calculate_model_image();
	}

	return logH;
}

double MyModel::log_likelihood() const
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
			if(sigma[i][j] < 1E100)
			{
				var = sigma[i][j]*sigma[i][j] + sigma0*sigma0
						+ sigma1*model_image[i][j];
				logL += -0.5*log(2.*M_PI*var) -
				0.5*pow(image[i][j] - model_image[i][j], 2)/var;
			}
		}
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	out<<setprecision(6);
	out<<' '<<sigma0<<' '<<sigma1<<' ';
	lens.print(out);
	source.print(out);

	// Make an image of the source (uses the ray resolution)
	const vector< vector<double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<double> >& y = Data::get_instance().get_y_rays();
	for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
			out<<source.evaluate(x[i][j], y[i][j], false)<<' ';

	for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
			out<<surface_brightness[i][j]<<' ';

	for(size_t i=0; i<model_image.size(); i++)
		for(size_t j=0; j<model_image[i].size(); j++)
			out<<model_image[i][j]<<' ';
}

string MyModel::description() const
{
	return string("lens parameters (b, q, rc, xc, yc, theta, shear, theta_shear), sigma0, sigma1, mock image");
}

void MyModel::shoot_rays(bool update)
{
	const vector< vector<double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<double> >& y = Data::get_instance().get_y_rays();

	double ax, ay;
	for(size_t i=0; i<xs.size(); i++)
	{
		for(size_t j=0; j<xs[i].size(); j++)
		{
			lens.alpha(x[i][j], y[i][j], ax, ay, update);
        
            if(update)
            {
                xs[i][j] -= ax;
                ys[i][j] -= ay;
            }
            else
            {
    			xs[i][j] = x[i][j] - ax;
        		ys[i][j] = y[i][j] - ay;
            }
		}
	}
}

void MyModel::calculate_surface_brightness(bool update)
{
    for(size_t i=0; i<xs.size(); i++)
    {
        for(size_t j=0; j<xs[i].size(); j++)
        {
            if(!update)
                surface_brightness[i][j] = 0.0;

            surface_brightness[i][j] += source.evaluate(xs[i][j], ys[i][j],
                                                                update);
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
}

void MyModel::calculate_model_image()
{
	int resolution = Data::get_instance().get_resolution();

	model_image.assign(Data::get_instance().get_ni(),
		vector<double>(Data::get_instance().get_nj(), 0.));

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

} // namespace Lensing2

