#ifndef Lensing2_Lens_h
#define Lensing2_Lens_h

#include <ostream>
#include "DNest4/code/RNG.h"

namespace Lensing2
{

class Lens
{
	protected:


	public:
		virtual ~Lens() { }

		// Deflection angle formula
		virtual void alpha(double x, double y,
					double& ax, double& ay) = 0;

		// MCMC related stuff
		virtual void from_prior(DNest4::RNG& rng) = 0;
		virtual double perturb(DNest4::RNG& rng) = 0;
		virtual void print(std::ostream& out) const = 0;
};

} // namespace Lensing2

#endif

