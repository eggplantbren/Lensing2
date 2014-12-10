#ifndef _MyModel3_
#define _MyModel3_

#include "Model.h"
#include "Sources/Sersic.h"
#include "Lenses/BlobbyNIE.h"
#include <vector>

namespace Lensing2
{

// FOR LENSED QSOS

class MyModel3:public DNest3::Model
{
	private:
		BlobbyNIE lens;

		// Numerically calculate the magnification at a given point
		double magnification(double x, double y) const;

	public:
		MyModel3();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

		static void test();
};

} // namespace Lensing2

#endif
