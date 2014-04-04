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
		// A maximum scale that we should use to judge a plausible range
		// for the parameters b and rc. Base this on the size of the
		// lensed image in whatever units are being used.
		const double scale;

		// Strength, axis ratio, core radius
		double b, q, rc;

		// Position
		double xc, yc;

		// Orientation angle
		double theta;

	public:
		PIEP(double scale);		

		// Needed methods
		void alpha(double x, double y, double& ax, double& ay);
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

}

#endif

