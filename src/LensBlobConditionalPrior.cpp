#include "LensBlobConditionalPrior.h"
#include "DNest4/code/RNG.h"
#include "DNest4/code/Utils.h"
#include <cmath>

using namespace DNest4;

const boost::math::normal LensBlobConditionalPrior::normal(0.0, 1.0);

LensBlobConditionalPrior::LensBlobConditionalPrior(double x_min, double x_max,
					double y_min, double y_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,size(sqrt((x_max - x_min)*(y_max - y_min)))
,xc(0.5*(x_min + x_max))
,yc(0.5*(y_min + y_max))
{

}

void LensBlobConditionalPrior::from_prior(RNG& rng)
{
    // mass \propto (einstein radius)^2
    // Don't want ER > size
    // =>
    // Don't want mass > size^2
	mu = exp(log(1E-6*size*size) + log(1E6)*rng.rand());

	b = exp(log(0.001*size) + log(300.0)*rng.rand());
	k = rng.rand();
	a = k*b;
}

double LensBlobConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;
	int which = rng.rand_int(3);

	if(which == 0)
	{
		mu = log(mu);
		mu += log(1E6)*rng.randh();
        DNest4::wrap(mu, log(1E-6*size*size), log(size*size));
		mu = exp(mu);
	}
	else if(which == 1)
	{
		b = log(b);
		b += log(300.0)*rng.randh();
        wrap(b, log(0.001*size), log(0.3*size));
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

double LensBlobConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    if(vec[2] < 0. || vec[3] < a || vec[3] > b)
        return -1E300;

    double logp = 0.0;

    double rsq = pow(vec[0] - xc, 2) + pow(vec[1] - yc, 2);
    logp += -0.5*rsq / (size * size);

    logp += -log(mu) - vec[2]/mu;
    logp += -log(b - a);

    return logp;
}

void LensBlobConditionalPrior::from_uniform(std::vector<double>& vec) const
{
	vec[0] = xc + size * quantile(normal, vec[0]);
	vec[1] = yc + size * quantile(normal, vec[1]);
	vec[2] = -mu*log(1. - vec[2]);
	vec[3] = a + (b - a)*vec[3];
}

void LensBlobConditionalPrior::to_uniform(std::vector<double>& vec) const
{
	vec[0] = cdf(normal, (vec[0] - xc) / size);
	vec[1] = cdf(normal, (vec[1] - yc) / size);
	vec[2] = 1. - exp(-vec[2]/mu);
	vec[3] = (vec[3] - a)/(b - a);
}

void LensBlobConditionalPrior::print(std::ostream& out) const
{
	out<<mu<<' '<<a<<' '<<b<<' ';
}

void LensBlobConditionalPrior::read(std::istream& in)
{
    in>>mu>>a>>b;
    k = a/b;
}

