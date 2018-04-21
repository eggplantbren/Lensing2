#include "BasicCircular.h"
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "DNest4/code/Distributions/Cauchy.h"
#include "DNest4/code/RNG.h"
#include "DNest4/code/Utils.h"
#include <cmath>

using namespace boost::math;
using namespace DNest4;

BasicCircular::BasicCircular(double x_min, double x_max,
					double y_min, double y_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,size(sqrt((x_max - x_min)*(y_max - y_min)))
{

}

void BasicCircular::from_prior(RNG& rng)
{
	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(rng.rand() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(rng.rand() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	width = exp(log(1E-2*size) + log(1E3)*rng.rand());
    shape = 0.1 + 1.9*rng.rand();

    DNest4::Cauchy cauchy(0.0, 5.0);
    do
    {
    	mu = cauchy.generate(rng);
    }while(std::abs(mu) > 50.0);
    mu = exp(mu);

	b = exp(log(1E-3*size) + log(1E3)*rng.rand());
	k = rng.rand();
	a = k*b;
}

double BasicCircular::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;
	int which = rng.rand_int(6);

	if(which == 0)
	{
		logH -= -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH -= -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));

		xc += (x_max - x_min)*rng.randh();
		yc += (y_max - y_min)*rng.randh();

		xc = mod(xc - x_min, x_max - x_min) + x_min;
		yc = mod(yc - y_min, y_max - y_min) + y_min;

		logH += -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH += -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));
	}
    else if(which == 1)
    {
        shape += 1.9*rng.randh();
        DNest4::wrap(shape, 0.1, 2.0);
    }
	else if(which == 2)
	{
		width = log(width);
		width += log(1E3)*pow(10., 1.5 - 6.*rng.rand())*rng.randn();
		width = mod(width - log(1E-2*size), log(1E3)) + log(1E-2*size);
		width = exp(width);
	}
	else if(which == 3)
	{
        DNest4::Cauchy cauchy(0.0, 5.0);

        mu = log(mu);
        logH += cauchy.perturb(mu, rng);
        if(std::abs(mu) > 50.0)
            return -1E300;
        mu = exp(mu);
	}
	else if(which == 4)
	{
		b = log(b);
		b += log(1E3)*pow(10., 1.5 - 6.*rng.rand())*rng.randn();
		b = mod(b - log(1E-3*size), log(1E3)) + log(1E-3*size);
		b = exp(b);
		a = k*b;
	}
	else
	{
		k += rng.randh();
		wrap(k, 0., 1.);
		a = k*b;
	}

	return logH;
}

double BasicCircular::log_pdf(const std::vector<double>& vec) const
{
    gamma_distribution<double> my_gamma(1.0/(shape*shape), width*shape*shape);

	if(vec[2] < 0. || vec[3] < a || vec[3] > b)
		return -1E300;

	double logp = 0.;
	double r = sqrt(pow(vec[0] - xc, 2) + pow(vec[1] - yc, 2));

    double alpha  = 1.0/(shape*shape);
    double lambda = 1.0/(width*shape*shape);

	logp += -log(r) - std::lgamma(alpha) + alpha*log(lambda)
                        + (alpha - 1.0)*log(r) - lambda*r;
	logp += -log(mu) - vec[2]/mu;
	logp += -log(b - a);

	return logp;
}

void BasicCircular::from_uniform(std::vector<double>& vec) const
{
    gamma_distribution<double> my_gamma(1.0/(shape*shape), width*shape*shape);

    double r;
    try
    {
        r = quantile(my_gamma, vec[0]);
    }
    catch(...)
    {
        return;
    }
	double phi = 2.*M_PI*vec[1];

	vec[0] = xc + r*cos(phi);
	vec[1] = yc + r*sin(phi);
	vec[2] = -mu*log(1. - vec[2]);
	vec[3] = a + (b - a)*vec[3];
}

void BasicCircular::to_uniform(std::vector<double>& vec) const
{
    gamma_distribution<double> my_gamma(1.0/(shape*shape), width*shape*shape);

	double r = sqrt(pow(vec[0] - xc, 2) + pow(vec[1] - yc, 2));
	double phi = atan2(vec[1] - yc, vec[0] - xc);
	if(phi < 0.)
		phi += 2.*M_PI;

    try
    {
      	vec[0] = cdf(my_gamma, r);
    }
    catch(...)
    {
        return;
    }

	vec[1] = phi/(2.*M_PI);
	vec[2] = 1. - exp(-vec[2]/mu);
	vec[3] = (vec[3] - a)/(b - a);
}

void BasicCircular::print(std::ostream& out) const
{
	out<<xc<<' '<<yc<<' '<<width<<' '<<shape<<' '<<mu<<' '<<a<<' '<<b<<' ';
}

void BasicCircular::read(std::istream& in)
{
    in>>xc>>yc>>width>>shape>>mu>>a>>b;
    k = a/b;
}


