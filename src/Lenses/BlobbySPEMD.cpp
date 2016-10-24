#include "BlobbySPEMD.h"

#include <cmath>
#include <cassert>

using namespace std;
using namespace DNest4;
using namespace Lensing2;

extern "C"
{
	void fastelldefl_(double* x, double* y, double* b, double* gam, double* q,
						double* rcsq, double alpha[]);
}

const bool BlobbySPEMD::disable_blobs = false;
const bool BlobbySPEMD::singular = true;

BlobbySPEMD::BlobbySPEMD(double x_min, double x_max, double y_min, double y_max)
:x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
,scale(sqrt((x_max - x_min)*(y_max - y_min)))
,blobs(4, 50, false,
	BasicUniform(x_min, x_max, y_min, y_max), PriorType::log_uniform)
,blobs_flag(false)
{
    if(x_max < x_min || y_max < y_min)
        throw std::logic_error("Invalid input to BlobbyNIE constructor.");
}

void BlobbySPEMD::alpha(double x, double y, double& ax, double& ay, bool update)
{
    ax = 0.0;
    ay = 0.0;

    if(!update)
    {
        // Rotate and center
        double xx =  (x - xc)*cos_theta + (y - yc)*sin_theta;
        double yy = -(x - xc)*sin_theta + (y - yc)*cos_theta;

        double aa[2];
        double rcsq = rc*rc;
        double coeff = 0.5*pow(bb, 2.*slope)*(2. - 2.*slope);

        fastelldefl_(&xx, &yy, &coeff, &slope, &qq, &rcsq, aa);
        double alphax = aa[0];
        double alphay = aa[1];

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
    }

	// Add blobs
	if(BlobbySPEMD::disable_blobs)
		return;

	const vector< vector<double> >& components = (update)?(blobs.get_added()):(blobs.get_components());
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

void BlobbySPEMD::from_prior(RNG& rng)
{
    while(true)
    {
	    b = rng.rand();
        if(rng.rand() <= (b/scale))
            break;
    }
	q = 0.05 + 0.95*rng.rand();

	// Stuff derived from b and q
	qq = q;
	if(qq == 1.)
		qq = 0.99999;
	bb = b/sqrt(qq); // Semi-major axis

	if(singular)
		rc = 1E-7*scale;
	else
		rc = exp(log(1E-3) + log(1E3)*rng.rand())*scale;

    slope = 0.1 + 0.8*rng.rand();

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

double BlobbySPEMD::perturb(RNG& rng)
{
	double logH = 0.;

	blobs_flag = false;
	if(rng.rand() <= 0.5 && !BlobbySPEMD::disable_blobs)
	{
		logH += blobs.perturb(rng);
		blobs_flag = true;
		return logH;
	}

	int which;
	do
	{
		which = rng.rand_int(8);
	}while(singular && which == 2);

	if(which == 0)
	{
        logH -= log(b/scale);
		b += rng.randh();
        DNest4::wrap(b, 0.0, scale);
        logH += log(b/scale);

		// Stuff derived from b and q
		qq = q;
		if(qq == 1.)
			qq = 0.99999;
		bb = b/sqrt(qq); // Semi-major axis
	}
	else if(which == 1)
	{
		q += 0.95*rng.randh();
		wrap(q, 0.05, 1.);

		// Stuff derived from b and q
		qq = q;
		if(qq == 1.)
			qq = 0.99999;
		bb = b/sqrt(qq); // Semi-major axis
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
		slope += 0.8*rng.randh();
		wrap(slope, 0.1, 0.9);
	}
	else if(which == 4)
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
	else if(which == 5)
	{
		theta += M_PI*rng.randh();
		theta = mod(theta, M_PI);
		cos_theta = cos(theta); sin_theta = sin(theta);
	}
	else if(which == 6)
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

void BlobbySPEMD::print(ostream& out) const
{
	out<<b<<' '<<q<<' '<<rc<<' '<<slope<<' '<<xc<<' '<<yc<<' '<<theta<<' ';
	out<<shear<<' '<<theta_shear<<' ';
	blobs.print(out); out<<' ';
}

void BlobbySPEMD::read(std::istream& in)
{
    in>>b>>q>>rc>>slope>>xc>>yc>>theta>>shear>>theta_shear;
    blobs.read(in);
}

