#ifndef Lensing2_Lens_h
#define Lensing2_Lens_h

#include <ostream>

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
		virtual void from_prior() = 0;
		virtual double perturb() = 0;
		virtual void print(std::ostream& out) const = 0;
};

} // namespace Lensing2

#endif

