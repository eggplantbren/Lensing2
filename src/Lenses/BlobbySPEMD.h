#ifndef Lensing2_BlobbySPEMD_h
#define Lensing2_BlobbySPEMD_h

#include <ostream>
#include <exception>
#include "../BasicUniform.h"
#include "../Lens.h"
#include "DNest4/code/RJObject/RJObject.h"

namespace Lensing2
{

/*
* A class for a single BlobbySPEMD.
*/
class BlobbySPEMD:public Lens
{
	private:
		static const bool singular;

		// Image dimensions
		double x_min, x_max, y_min, y_max, scale;

		// Strength, axis ratio, core radius
		double b, q, rc, slope;

		// Derived parameters based on b and q
		double bb, qq;

		// Position
		double xc, yc;

		// Orientation angle
		double theta, cos_theta, sin_theta;

		// External shear
		double shear;
		double theta_shear, cos_theta_shear, sin_theta_shear;

		DNest4::RJObject<BasicUniform> blobs;
		bool blobs_flag;

		// A flag to disable the blobs
		static const bool disable_blobs;

	public:
		// Constructor: pass in dimensions of image
		BlobbySPEMD(double x_min, double x_max, double y_min, double y_max);

		// Needed methods
		void alpha(double x, double y, double& ax, double& ay, bool update=false);
		void from_prior(DNest4::RNG& rng);
		double perturb(DNest4::RNG& rng);
		void print(std::ostream& out) const;

		// Getter
		const DNest4::RJObject<BasicUniform>& get_blobs() const
        { return blobs; }

		// If true, the most recent perturb only involved the blobs
		bool get_blobs_flag() const { return blobs_flag; }
};

}

#endif

