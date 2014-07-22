#include "BlobbyNIE.h"

#include "RandomNumberGenerator.h"
#include "Utils.h"

#include <cmath>
#include <cassert>

using namespace std;
using namespace DNest3;
using namespace Lensing2;

BlobbyNIE::BlobbyNIE(double x_min, double x_max, double y_min, double y_max)
:x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
,scale(sqrt((x_max - x_min)*(y_max - y_min)))
,blobs(4, 10, false,
	BasicCircular(x_min, x_max, y_min, y_max, 1E-4*scale, 10*scale))
,blobs_flag(false)
{
	assert(x_max > x_min && y_max > y_min);
}

void BlobbyNIE::alpha(double x, double y, double& ax, double& ay) const
{
	// Rotate and center
	double xx =  (x - xc)*cos_theta + (y - yc)*sin_theta;
	double yy = -(x - xc)*sin_theta + (y - yc)*cos_theta;

	double qq = q;

	// Based on Keeton and Kochanek (1998) equations 6 and 7.
	// Might be better to symmetrise but I haven't done it yet.
	if(qq == 1.)
		qq = 0.99999;

	// Formulae use the minor axis length
	double bb = b*sqrt(qq);

	double psi = sqrt(qq*qq*(xx*xx + rc*rc) + yy*yy);
	double q_term = sqrt(1. - qq*qq);
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

void BlobbyNIE::alpha_diff(double x, double y, double& ax, double& ay) const
{
	ax = 0.;
	ay = 0.;

	// Add blobs
	const vector< vector<double> >& added = blobs.get_added();
	double rsq, widthsq, Menc;
	for(size_t i=0; i<added.size(); i++)
	{
		rsq = pow(x - added[i][0], 2)
				+ pow(y - added[i][1], 2);
		widthsq = pow(added[i][3], 2);

		if(rsq < widthsq)
		{
			Menc = 4.*added[i][2]/widthsq*(0.5*rsq -
				rsq*rsq/(4*widthsq));
		}
		else
			Menc = added[i][2];
		ax += Menc*(x - added[i][0])/(M_PI*rsq);
		ay += Menc*(y - added[i][1])/(M_PI*rsq);
	}

	// Remove blobs
	const vector< vector<double> >& removed = blobs.get_removed();
	for(size_t i=0; i<removed.size(); i++)
	{
		rsq = pow(x - removed[i][0], 2)
				+ pow(y - removed[i][1], 2);
		widthsq = pow(removed[i][3], 2);

		if(rsq < widthsq)
		{
			Menc = 4.*removed[i][2]/widthsq*(0.5*rsq -
				rsq*rsq/(4*widthsq));
		}
		else
			Menc = removed[i][2];
		ax -= Menc*(x - removed[i][0])/(M_PI*rsq);
		ay -= Menc*(y - removed[i][1])/(M_PI*rsq);
	}
}


void BlobbyNIE::from_prior()
{
	b = exp(log(1E-3) + log(1E3)*randomU())*scale;
	q = exp(log(0.1) + log(10.)*randomU());
	rc = exp(log(1E-3) + log(1E3)*randomU())*scale;

	do
	{
		xc = 0.5*(x_max + x_min) +
			0.1*(x_max - x_min)*tan(M_PI*(randomU() - 0.5));
		yc = 0.5*(y_max + y_min) +
			0.1*(y_max - y_min)*tan(M_PI*(randomU() - 0.5));
	}while(xc < x_min || xc > x_max || yc < y_min || yc > y_max);

	theta = M_PI*randomU();
	cos_theta = cos(theta); sin_theta = sin(theta);

	// Half-cauchy prior
	shear = 0.05*tan(M_PI*(0.5*randomU()));
	theta_shear = 2.*M_PI*randomU();
	cos_theta_shear = cos(theta_shear); sin_theta_shear = sin(theta_shear);

	blobs.fromPrior();
}

double BlobbyNIE::perturb()
{
	double logH = 0.;

	blobs_flag = false;
	if(randomU() <= 0.5)
	{
		logH += blobs.perturb();
		blobs_flag = true;
		return logH;
	}

	int which = randInt(7);

	if(which == 0)
	{
		b = log(b/scale);
		b += log(1E3)*randh();
		b = mod(b - log(1E-3), log(1E3)) + log(1E-3);
		b = scale*exp(b);
	}
	else if(which == 1)
	{
		q = log(q);
		q += log(10.)*randh();
		q = mod(q - log(0.1), log(10.)) + log(0.1);
		q = exp(q);
	}
	else if(which == 2)
	{
		rc = log(rc/scale);
		rc += log(1E3)*randh();
		rc = mod(rc - log(1E-3), log(1E3)) + log(1E-3);
		rc = scale*exp(rc);
	}
	else if(which == 3)
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
	else if(which == 4)
	{
		theta += M_PI*randh();
		theta = mod(theta, M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}
	else if(which == 5)
	{
		shear = atan(shear/0.05)/M_PI/0.5;
		shear += randh();
		shear = mod(shear, 1.);
		shear = 0.05*tan(M_PI*(0.5*shear));
	}
	else
	{
		theta_shear += 2.*M_PI*randh();
		theta_shear = mod(theta_shear, 2.*M_PI);
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

