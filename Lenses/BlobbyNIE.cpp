#include "BlobbyNIE.h"

#include <cmath>
#include <cassert>

using namespace std;
using namespace DNest4;
using namespace Lensing2;

const bool BlobbyNIE::disable_blobs = false;
const bool BlobbyNIE::singular = true;

BlobbyNIE::BlobbyNIE(double x_min, double x_max, double y_min, double y_max)
:x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
,scale(sqrt((x_max - x_min)*(y_max - y_min)))
,blobs(4, 10, false,
	BasicCircular(x_min, x_max, y_min, y_max), PriorType::log_uniform)
,blobs_flag(false)
{
	if(x_max < x_min || y_max < y_min)
        throw std::logic_error("Invalid input to BlobbyNIE constructor.");
}

void BlobbyNIE::alpha(double x, double y, double& ax, double& ay)
{
	// Rotate and center
	double xx =  (x - xc)*cos_theta + (y - yc)*sin_theta;
	double yy = -(x - xc)*sin_theta + (y - yc)*cos_theta;

	double psi = sqrt(qq*qq*(xx*xx + rc*rc) + yy*yy);
	double alphax = bb/q_term*atan(q_term*xx/(psi + rc));
	double alphay = bb/q_term*atanh(q_term*yy/(psi + qq*qq*rc));

	// Rotate back
	ax = alphax*cos_theta - alphay*sin_theta;
	ay = alphax*sin_theta + alphay*cos_theta;

	// Go into shear coordinate system
	xx =  x*cos_theta_shear + y*sin_theta_shear;
	yy = -x*sin_theta_shear + y*cos_theta_shear;

	// Calculate external shear
	alphax = -shear*xx;
	alphay = shear*yy;

	// Add external shear
	ax += alphax*cos_theta_shear - alphay*sin_theta_shear;
	ay += alphax*sin_theta_shear + alphay*cos_theta_shear;

	// Add blobs
	if(BlobbyNIE::disable_blobs)
		return;

	const vector< vector<double> >& components = blobs.get_components();
	double rsq, widthsq, Menc;
	for(size_t i=0; i<components.size(); i++)
	{
		rsq = pow(x - components[i][0], 2)
				+ pow(y - components[i][1], 2);
		widthsq = pow(components[i][3], 2);

		if(rsq < widthsq)
		{
			Menc = 4.*components[i][2]/widthsq*(0.5*rsq -
				rsq*rsq/(4*widthsq));
		}
		else
			Menc = components[i][2];
		ax += Menc*(x - components[i][0])/(M_PI*rsq);
		ay += Menc*(y - components[i][1])/(M_PI*rsq);
	}
}

void BlobbyNIE::from_prior(RNG& rng)
{
	b = exp(log(1E-3) + log(1E3)*rng.rand())*scale;
	q = 0.05 + 0.95*rng.rand();

	// Stuff derived from b and q
	qq = q;
	if(qq == 1.)
		qq = 0.99999;
	q_term = sqrt(1. - qq*qq);
	bb = b*sqrt(qq); // Minor axis

	if(singular)
		rc = 1E-7*scale;
	else
		rc = exp(log(1E-3) + log(1E3)*rng.rand())*scale;

	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(rng.rand() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(rng.rand() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	theta = M_PI*rng.rand();
	cos_theta = cos(theta); sin_theta = sin(theta);

	// Half-cauchy prior
	shear = 0.05*tan(M_PI*(0.5*rng.rand()));
	theta_shear = M_PI*rng.rand();
	cos_theta_shear = cos(theta_shear); sin_theta_shear = sin(theta_shear);

	blobs.from_prior(rng);
}

double BlobbyNIE::perturb(RNG& rng)
{
	double logH = 0.;

	blobs_flag = false;
	if(rng.rand() <= 0.5 && !BlobbyNIE::disable_blobs)
	{
		logH += blobs.perturb(rng);
		blobs_flag = true;
		return logH;
	}

	int which;
	do
	{
		which = rng.rand_int(7);
	}while(singular && which == 2);

	if(which == 0)
	{
		b = log(b/scale);
		b += log(1E3)*rng.randh();
		b = mod(b - log(1E-3), log(1E3)) + log(1E-3);
		b = scale*exp(b);

		// Stuff derived from b and q
		qq = q;
		if(qq == 1.)
			qq = 0.99999;
		q_term = sqrt(1. - qq*qq);
		bb = b*sqrt(qq); // Minor axis
	}
	else if(which == 1)
	{
		q += 0.95*rng.randh();
		wrap(q, 0.05, 1.);

		// Stuff derived from b and q
		qq = q;
		if(qq == 1.)
			qq = 0.99999;
		q_term = sqrt(1. - qq*qq);
		bb = b*sqrt(qq); // Minor axis
	}
	else if(which == 2)
	{
		rc = log(rc/scale);
		rc += log(1E3)*rng.randh();
		rc = mod(rc - log(1E-3), log(1E3)) + log(1E-3);
		rc = scale*exp(rc);
	}
	else if(which == 3)
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
	else if(which == 4)
	{
		theta += M_PI*rng.randh();
		theta = mod(theta, M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}
	else if(which == 5)
	{
		shear = atan(shear/0.05)/M_PI/0.5;
		shear += rng.randh();
		shear = mod(shear, 1.);
		shear = 0.05*tan(M_PI*(0.5*shear));
	}
	else
	{
		theta_shear += M_PI*rng.randh();
		theta_shear = mod(theta_shear, M_PI);
		cos_theta_shear = cos(theta_shear); sin_theta_shear = sin(theta_shear);
	}

	return logH;
}

void BlobbyNIE::print(ostream& out) const
{
	out<<b<<' '<<q<<' '<<rc<<' '<<xc<<' '<<yc<<' '<<theta<<' ';
	out<<shear<<' '<<theta_shear<<' ';
	blobs.print(out); out<<' ';
}

int BlobbyNIE::get_size_of_diff() const
{
	return static_cast<int>(blobs.get_added().size()
					+ blobs.get_removed().size());
}

int BlobbyNIE::get_num_components() const
{
	return static_cast<int>(blobs.get_components().size());
}

