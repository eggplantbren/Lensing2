#include "BasicCircular.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace DNest3;

BasicCircular::BasicCircular(double x_min, double x_max,
					double y_min, double y_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,size(sqrt((x_max - x_min)*(y_max - y_min)))
{

}

void BasicCircular::fromPrior()
{
	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(randomU() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(randomU() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	width = exp(log(1E-2*size) + log(1E3)*randomU());

	do
	{
		mu = tan(M_PI*(randomU() - 0.5));
	}while(abs(mu) > 20.);
	mu = exp(mu);

	mu_w = exp(log(1E-3*size) + log(1E3)*randomU());
}

double BasicCircular::perturb_parameters()
{
	double logH = 0.;
	int which = randInt(4);

	if(which == 0)
	{
		logH -= -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH -= -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));

		xc += (x_max - x_min)*randh();
		yc += (y_max - y_min)*randh();

		xc = mod(xc - x_min, x_max - x_min) + x_min;
		yc = mod(yc - y_min, y_max - y_min) + y_min;

		logH += -log(1. + pow((xc - 0.5*(x_min + x_max))/(0.1*(x_max - x_min)), 2));
		logH += -log(1. + pow((yc - 0.5*(y_min + y_max))/(0.1*(y_max - y_min)), 2));
	}
	else if(which == 1)
	{
		width = log(width);
		width += log(1E3)*pow(10., 1.5 - 6.*randomU())*randn();
		width = mod(width - log(1E-2*size), log(1E3)) + log(1E-2*size);
		width = exp(width);
	}
	else if(which == 2)
	{
		mu = log(mu);
		mu = atan(mu)/M_PI + 0.5;
		mu += randh();
		mu = mod(mu, 1.);
		mu = tan(M_PI*(mu - 0.5));
		if(abs(mu) > 20.)
			logH = -1E300;
		mu = exp(mu);
	}
	else
	{
		mu_w = log(mu_w);
		mu_w += log(1E3)*pow(10., 1.5 - 6.*randomU())*randn();
		mu_w = mod(mu_w - log(1E-3*size), log(1E3)) + log(1E-3*size);
		mu_w = exp(mu_w);
	}

	return logH;
}

double BasicCircular::log_pdf(const std::vector<double>& vec) const
{
	if(vec[2] < 0. || vec[3] < 0.)
		return -1E300;

	double logp = 0.;
	double r = sqrt(pow(vec[0] - xc, 2) + pow(vec[1] - yc, 2));

	logp += -log(r) - log(width) - r/width;
	logp += -log(mu) - vec[2]/mu;
	logp += -log(mu_w) - vec[3]/mu_w;

	return logp;
}

void BasicCircular::from_uniform(std::vector<double>& vec) const
{
	double r = -width*log(1. - vec[0]);
	double phi = 2.*M_PI*vec[1];

	vec[0] = xc + r*cos(phi);
	vec[1] = yc + r*sin(phi);
	vec[2] = -mu*log(1. - vec[2]);
	vec[3] = -mu_w*log(1. - vec[3]);
}

void BasicCircular::to_uniform(std::vector<double>& vec) const
{
	double r = sqrt(pow(vec[0] - xc, 2) + pow(vec[1] - yc, 2));
	double phi = atan2(vec[1] - yc, vec[0] - xc);
	if(phi < 0.)
		phi += 2.*M_PI;

	vec[0] = 1. - exp(-r/width);
	vec[1] = phi/(2.*M_PI);
	vec[2] = 1. - exp(-vec[2]/mu);
	vec[3] = 1. - exp(-vec[3]/mu_w);
}

void BasicCircular::print(std::ostream& out) const
{
	out<<xc<<' '<<yc<<' '<<width<<' '<<mu<<' '<<mu_w<<' ';
}

