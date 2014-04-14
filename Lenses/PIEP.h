#ifndef Lensing2_PIEP_h
#define Lensing2_PIEP_h

#include <ostream>
#include "Lens.h"

namespace Lensing2
{

/*
* A class for a single PIEP.
*/
class PIEP:public Lens
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

	public:
		// Constructor: pass in dimensions of image
		PIEP(double x_min, double x_max, double y_min, double y_max);

		// Needed methods
		void alpha(double x, double y, double& ax, double& ay) const;
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

}

#endif

