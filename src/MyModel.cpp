#include "MyModel.h"
#include "Data.h"
#include "DNest4/code/Distributions/Cauchy.h"
#include <cmath>
#include <sstream>

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
    psf_power = 2*rng.rand();


    DNest4::Cauchy cauchy(0.0, 5.0);

    do
    {
    	sigma0 = cauchy.generate(rng);
    }while(std::abs(sigma0) > 50.0);
    sigma0 = exp(sigma0);

    do
    {
    	sigma1 = cauchy.generate(rng);
    }while(std::abs(sigma1) > 50.0);
    sigma1 = exp(sigma1);

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
        int which = rng.rand_int(3);

        if(which == 0)
        {
            DNest4::Cauchy cauchy(0.0, 5.0);

            sigma0 = log(sigma0);
            logH += cauchy.perturb(sigma0, rng);
            if(std::abs(sigma0) > 50.0)
                return -1E300;
            sigma0 = exp(sigma0);
        }
        else if(which == 1)
        {
            DNest4::Cauchy cauchy(0.0, 5.0);

            sigma1 = log(sigma1);
            logH += cauchy.perturb(sigma1, rng);
            if(std::abs(sigma1) > 50.0)
                return -1E300;
            sigma1 = exp(sigma1);
        }
        else
        {
            psf_power += 2*rng.rand();
            wrap(psf_power, 0.0, 2.0);
            if(Data::get_instance().psf_is_highres())
	            calculate_surface_brightness();
	        calculate_model_image();
        }
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
	out<<' '<<sigma0<<' '<<sigma1<<' '<<psf_power<<' ';
	lens.print(out);
	source.print(out);

	// Make an image of the source (uses the ray resolution)
	const vector< vector<double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<double> >& y = Data::get_instance().get_y_rays();
    std::vector<std::vector<double>> evals
            (x.size(), std::vector<double>(x[0].size()));
    source.evaluate(x, y, evals, false);

	for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
			out<<evals[i][j]<<' ';

	for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
			out<<surface_brightness[i][j]<<' ';

	for(size_t i=0; i<model_image.size(); i++)
		for(size_t j=0; j<model_image[i].size(); j++)
			out<<model_image[i][j]<<' ';
}

string MyModel::description() const
{
    stringstream s;
    s<<"sigma0, sigma1, psf_power, ";
    s<<"b, q, rc, slope, xc, yc, theta, shear, theta_shear, ";
    s<<"dim_lens_blobs, max_num_lens_blobs, ";
    s<<"mu_lens_blobs, a_lens_blobs, b_lens_blobs, ";
    s<<"num_lens_blobs, ";

    for(int i=0; i<50; ++i)
        s<<"lens_blob_x["<<i<<"], ";
    for(int i=0; i<50; ++i)
        s<<"lens_blob_y["<<i<<"], ";
    for(int i=0; i<50; ++i)
        s<<"lens_blob_mass["<<i<<"], ";
    for(int i=0; i<50; ++i)
        s<<"lens_blob_width["<<i<<"], ";

    s<<"dim_source_blobs, max_num_source_blobs, ";
    s<<"xc_source_blobs, yc_source_blobs, width_source_blobs, ";
    s<<"mu_source_blobs, a_source_blobs, b_source_blobs, ";
    s<<"num_source_blobs, ";

    for(int i=0; i<500; ++i)
        s<<"source_blob_x["<<i<<"], ";
    for(int i=0; i<500; ++i)
        s<<"source_blob_y["<<i<<"], ";
    for(int i=0; i<500; ++i)
        s<<"source_blob_mass["<<i<<"], ";
    for(int i=0; i<500; ++i)
        s<<"source_blob_width["<<i<<"], ";

    for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
            s<<"source["<<i<<"]["<<j<<"], ";

    for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
            s<<"image["<<i<<"]["<<j<<"], ";

    for(size_t i=0; i<xs.size(); i++)
		for(size_t j=0; j<xs[i].size(); j++)
            s<<"blurred_image["<<i<<"]["<<j<<"], ";

	return s.str();
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
    std::vector<std::vector<double>> delta
            (xs.size(), std::vector<double>(xs[0].size()));

    source.evaluate(xs, ys, delta, update);

    if(Data::get_instance().psf_is_highres())
    {
        // Blur using the PSF
        const PSF& psf = Data::get_instance().get_psf();
        auto psf2 = psf;
        psf2.calculate_fft(delta.size(), delta[0].size(), psf_power);
        psf2.blur_image2(delta);
    }

    // Increment or replace surface_brightness
    if(update)
    {
        for(size_t i=0; i<xs.size(); ++i)
            for(size_t j=0; j<xs[i].size(); ++j)
                surface_brightness[i][j] += delta[i][j];  
    }
    else
        surface_brightness = delta;
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
        auto psf2 = psf;
        psf2.calculate_fft(model_image.size(),
                            model_image.size(), psf_power);
		psf2.blur_image2(model_image);
	}
}

void MyModel::read(std::istream& in)
{
    in>>sigma0>>sigma1>>psf_power;
    lens.read(in);
    source.read(in);

	shoot_rays();
	calculate_surface_brightness();
	calculate_model_image();
}

} // namespace Lensing2

