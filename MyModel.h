#ifndef _MyModel_
#define _MyModel_

#include "Model.h"
#include "Sources/Blobby.h"
#include "Lenses/PIEP.h"
#include <vector>

namespace Lensing2
{

class MyModel:public DNest3::Model
{
	private:
		Blobby source;
		PIEP lens;

	public:
		MyModel();

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
};

} // namespace Lensing2

#endif

