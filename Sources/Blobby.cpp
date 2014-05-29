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

	const vector< vector<double> >& components = blobs.get_components();

	double rsq, widthsq;
	for(size_t i=0; i<components.size(); i++)
	{
		rsq = pow(x - components[i][0], 2)
				+ pow(y - components[i][1], 2);
		widthsq = pow(components[i][3], 2);

		if(rsq < widthsq)
			f += components[i][2]*
				2./M_PI*(1. - rsq/widthsq)/widthsq;
	}

	return f;
}

int Blobby::get_size_of_diff() const
{
	return static_cast<int>(blobs.get_added().size()
					+ blobs.get_removed().size());
}

int Blobby::get_num_components() const
{
	return static_cast<int>(blobs.get_components().size());
}


double Blobby::evaluate_diff(double x, double y) const
{
	double f = 0.;

	// Added components
	const vector< vector<double> >& added = blobs.get_added();

	double rsq, widthsq;
	for(size_t i=0; i<added.size(); i++)
	{
		rsq = pow(x - added[i][0], 2)
				+ pow(y - added[i][1], 2);
		widthsq = pow(added[i][3], 2);

		if(rsq < widthsq)
			f += added[i][2]*
				2./M_PI*(1. - rsq/widthsq)/widthsq;
	}

	// Subtract removed components
	const vector< vector<double> >& removed = blobs.get_removed();
	for(size_t i=0; i<removed.size(); i++)
	{
		rsq = pow(x - removed[i][0], 2)
				+ pow(y - removed[i][1], 2);
		widthsq = pow(removed[i][3], 2);

		if(rsq < widthsq)
			f -= removed[i][2]*
				2./M_PI*(1. - rsq/widthsq)/widthsq;
	}

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
	blobs.print(out);
}

