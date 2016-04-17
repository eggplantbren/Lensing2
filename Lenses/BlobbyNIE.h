#ifndef Lensing2_BlobbyNIE_h
#define Lensing2_BlobbyNIE_h

#include <ostream>
#include <exception>
#include "../BasicCircular.h"
#include "../Lens.h"
#include "DNest4/code/RJObject/RJObject.h"

namespace Lensing2
{

/*
* A class for a single BlobbyNIE.
*/
class BlobbyNIE:public Lens
{
	private:
		static const bool singular;

		// Image dimensions
		double x_min, x_max, y_min, y_max, scale;

		// Strength, axis ratio, core radius
		double b, q, rc;

		// Derived parameters based on b and q
		double bb, qq, q_term;

		// Position
		double xc, yc;

		// Orientation angle
		double theta, cos_theta, sin_theta;

		// External shear
		double shear;
		double theta_shear, cos_theta_shear, sin_theta_shear;

		DNest4::RJObject<BasicCircular> blobs;
		bool blobs_flag;

		// A flag to disable the blobs
		static const bool disable_blobs;

	public:
		// Constructor: pass in dimensions of image
		BlobbyNIE(double x_min, double x_max, double y_min, double y_max);

		// Needed methods
		void alpha(double x, double y, double& ax, double& ay);
		void alpha_diff(double x, double y, double& ax, double& ay) const;
		void from_prior(DNest4::RNG& rng);
		double perturb(DNest4::RNG& rng);
		void print(std::ostream& out) const;

		// Methods for partial update
		int get_size_of_diff() const;
		int get_num_components() const;

		// If true, the most recent perturb only involved the blobs
		bool get_blobs_flag() const { return blobs_flag; }
};

}

#endif

