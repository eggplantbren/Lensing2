#ifndef Lensing2_Sersic_h
#define Lensing2_Sersic_h

#include "Source.h"
#include "../BasicCircular.h"
#include <RJObject.h>

namespace Lensing2
{

class Sersic:public Source
{
	private:
		// Boundaries of the image
		double x_min, x_max, y_min, y_max;
		// Image size
		double image_size;

		// Sersic parameters
		double Ie, m, Re;
		double one_over_Re_squared;
		double bm; // Related to m

		// Central position, axis ratio, angle
		double xc, yc, q, theta;

		// Stuff related to theta
		double cos_theta, sin_theta;

	public:
		// Pass in image scale and flux scale
		Sersic(double x_min, double x_max,
					double y_min, double y_max);

		// Required methods
		double evaluate(double x, double y) const;

		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

} // namespace Lensing2

#endif

