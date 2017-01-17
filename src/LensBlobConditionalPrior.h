#ifndef Lensing2_LensBlobConditionalPrior
#define Lensing2_LensBlobConditionalPrior

#include "DNest4/code/RJObject/ConditionalPriors/ConditionalPrior.h"
#include <istream>
#include <boost/math/distributions/normal.hpp>

class LensBlobConditionalPrior:public DNest4::ConditionalPrior
{
	private:
		// Limits
		double x_min, x_max, y_min, y_max, size;
        double xc, yc;

		// Mean of exponential distribution for masses
		double mu;

		// Uniform distribution for widths
		double b, k, a;

        // A standard normal distribution
        static const boost::math::normal normal;

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		LensBlobConditionalPrior(double x_min, double x_max,
					double y_min, double y_max);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
        void read(std::istream& in);
};

#endif

