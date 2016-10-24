#ifndef Lensing2_BasicUniform
#define Lensing2_BasicUniform

#include "DNest4/code/RJObject/ConditionalPriors/ConditionalPrior.h"

class BasicUniform:public DNest4::ConditionalPrior
{
	private:
		// Limits
		double x_min, x_max, y_min, y_max, size;

		// Mean of exponential distribution for masses
		double mu;

		// Uniform distribution for widths
		double b, k, a;

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		BasicUniform(double x_min, double x_max,
					double y_min, double y_max);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
};

#endif

