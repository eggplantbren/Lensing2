#include "Blobby.h"
#include "DNest4/code/RNG.h"

using namespace std;
using namespace DNest4;

namespace Lensing2
{

Blobby::Blobby(double x_min, double x_max,
					double y_min, double y_max)
:blobs(4, 500, false,
	BasicCircular(x_min, x_max, y_min, y_max), PriorType::log_uniform)
{

}

double Blobby::evaluate(double x, double y, bool update) const
{
	double f = 0.0;

	const vector< vector<double> >& components = (update)?
                                                 (blobs.get_added()):
                                                 (blobs.get_components());

	double rsq, widthsq, one_over;
    double c = 2.0/M_PI;

	for(size_t i=0; i<components.size(); ++i)
	{
		rsq = pow(x - components[i][0], 2) + pow(y - components[i][1], 2);
		widthsq = components[i][3]*components[i][3];
		one_over = 1.0/widthsq;

		if(rsq < widthsq)
			f += components[i][2]*c*(1.0 - rsq*one_over)*one_over;
	}

	return f;
}

void Blobby::from_prior(RNG& rng)
{
	blobs.from_prior(rng);
}

double Blobby::perturb(RNG& rng)
{
	double logH = 0.;

	logH += blobs.perturb(rng);

	return logH;
}

void Blobby::print(ostream& out) const
{
	blobs.print(out);
}

void Blobby::read(istream& in)
{
    blobs.read(in);
}

} // namespace Lensing2

