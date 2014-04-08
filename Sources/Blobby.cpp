#include "Blobby.h"
#include "RandomNumberGenerator.h"

using namespace std;
using namespace Lensing2;
using namespace DNest3;

Blobby::Blobby(double x_min, double x_max,
					double y_min, double y_max,
					double mu_min, double mu_max)
:blobs(4, 100, false,
	BasicCircular(x_min, x_max, y_min, y_max, mu_min, mu_max))
{

}

double Blobby::evaluate(double x, double y) const
{
	double f = 0.;

	f = x + y;

	return f;
}

void Blobby::from_prior()
{
	blobs.fromPrior();
}

double Blobby::perturb()
{
	double logH = 0.;

	logH += blobs.perturb();

	return logH;
}

void Blobby::print(ostream& out) const
{
	out<<' ';
}

