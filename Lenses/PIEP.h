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
		// Strength, axis ratio, core radius
		double b, q, rc;

		// Position
		double xc, yc;

		// Orientation angle
		double theta;

	public:
		// Needed methods
		void alpha(double x, double y, double& ax, double& ay);
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

}

#endif

