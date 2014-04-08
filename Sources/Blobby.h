#ifndef Lensing2_Blobby_h
#define Lensing2_Blobby_h

#include "Source.h"
#include "BasicCircular.h"
#include <RJObject.h>

namespace Lensing2
{

class Blobby:public Source
{
	private:
		RJObject<BasicCircular> blobs;

	public:
		// Pass in image scale and flux scale
		Blobby(double x_min, double x_max,
					double y_min, double y_max,
					double mu_min, double mu_max);

		// Required methods
		double evaluate(double x, double y) const;
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

} // namespace Lensing2

#endif

