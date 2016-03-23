#ifndef Lensing2_Source
#define Lensing2_Source

#include <ostream>
#include "DNest4/code/RNG.h"

namespace Lensing2
{

class Source
{
	protected:


	public:
		virtual ~Source() { }

		// Evaluate at a single position
		virtual double evaluate(double x, double y) const = 0;

		// MCMC related stuff
		virtual void from_prior(DNest4::RNG& rng) = 0;
		virtual double perturb(DNest4::RNG& rng) = 0;

		virtual void print(std::ostream& out) const = 0;
};

} // namespace Lensing2

#endif

