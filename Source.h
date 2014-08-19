#ifndef Lensing2_Source_h
#define Lensing2_Source_h

#include <ostream>

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
		virtual void from_prior() = 0;
		virtual double perturb() = 0;

		virtual void print(std::ostream& out) const = 0;
};

} // namespace Lensing2

#endif

