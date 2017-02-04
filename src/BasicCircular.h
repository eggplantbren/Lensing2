#ifndef Lensing2_BasicCircular
#define Lensing2_BasicCircular

#include "DNest4/code/RJObject/ConditionalPriors/ConditionalPrior.h"
#include <istream>

class BasicCircular:public DNest4::ConditionalPrior
{
	private:
		// Limits
		double x_min, x_max, y_min, y_max, size;

		// Center and width of circle
		double xc, yc;
		double width;
        double shape;

		// Mean of exponential distribution for masses
		double mu;

		// Uniform distribution for widths
		double b, k, a;

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		BasicCircular(double x_min, double x_max,
					  double y_min, double y_max);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
        void read(std::istream& in);
};

#endif

