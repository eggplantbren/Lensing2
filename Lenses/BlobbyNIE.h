#ifndef Lensing2_BlobbyNIE_h
#define Lensing2_BlobbyNIE_h

#include <ostream>
#include "../BasicCircular.h"
#include <RJObject.h>
#include "Lens.h"

namespace Lensing2
{

/*
* A class for a single BlobbyNIE.
*/
class BlobbyNIE:public Lens
{
	private:
		// Image dimensions
		double x_min, x_max, y_min, y_max, scale;

		// Strength, axis ratio, core radius
		double b, q, rc;

		// Position
		double xc, yc;

		// Orientation angle
		double theta, cos_theta, sin_theta;

		// External shear
		double shear;
		double theta_shear, cos_theta_shear, sin_theta_shear;

		RJObject<BasicCircular> blobs;

	public:
		// Constructor: pass in dimensions of image
		BlobbyNIE(double x_min, double x_max, double y_min, double y_max);

		// Needed methods
		void alpha(double x, double y, double& ax, double& ay) const;
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;

		// Methods for partial update
		int get_size_of_diff() const;
		int get_num_components() const;
};

}

#endif

