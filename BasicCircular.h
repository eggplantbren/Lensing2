#ifndef _BasicCircular_
#define _BasicCircular_

#include <Distributions/Distribution.h>

class BasicCircular:public Distribution
{
	private:
		// Limits
		double x_min, x_max, y_min, y_max, size;

		// Center and width of circle
		double xc, yc;
		double width;

		// Mean of exponential distribution for masses
		double mu;

		// Mean of exponential distribution for blob widths
		double mu_w;

		double perturb_parameters();

	public:
		BasicCircular(double x_min, double x_max,
					double y_min, double y_max);

		void fromPrior();

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
};

#endif

