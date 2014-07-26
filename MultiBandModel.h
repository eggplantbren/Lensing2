#ifndef _MultiBandModel_
#define _MultiBandModel_

#include "Model.h"
#include "Sources/Blobby.h"
#include "Lenses/BlobbyNIE.h"
#include <vector>

namespace Lensing2
{

class MultiBandModel:public DNest3::Model
{
	private:
		Blobby source;
		BlobbyNIE lens;

		// Noise model - constant variance plus component
		// that depends on the flux of the model
		double sigma0, sigma1;

		// Source plane position of rays
		std::vector< std::vector<long double> > xs;
		std::vector< std::vector<long double> > ys;

		// Surface brightness of rays
		std::vector< std::vector<long double> > surface_brightness;

		// Model image
		std::vector< std::vector<long double> > model_image;

		void shoot_rays();
		void update_rays();
		void calculate_surface_brightness();
		void update_surface_brightness();
		void calculate_model_image();

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

