#ifndef _MultiBandModel_
#define _MultiBandModel_

#include "Model.h"
#include "MyModel.h"
#include <vector>

namespace Lensing2
{

class MultiBandModel:public DNest3::Model
{
	private:
		std::vector<MyModel> bands;

	public:
		MultiBandModel();

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

