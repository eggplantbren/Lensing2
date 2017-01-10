#ifndef Lensing2_MyModel
#define Lensing2_MyModel

#include "Sources/Blobby.h"
#include "Lenses/BlobbySPEMD.h"
#include <vector>
#include <istream>

namespace Lensing2
{

class MyModel
{
	private:
		Blobby source;
		BlobbySPEMD lens;

        double psf_power;
        std::vector<std::vector<double>> alt_psf;

		// Noise model - constant variance plus component
		// that depends on the flux of the model
		double sigma0, sigma1;

		// Source plane position of rays
		std::vector< std::vector<double> > xs;
		std::vector< std::vector<double> > ys;

		// Surface brightness of rays
		std::vector< std::vector<double> > surface_brightness;

		// Model image
		std::vector< std::vector<double> > model_image;

		void shoot_rays(bool update=false);
		void calculate_surface_brightness(bool update=false);
		void calculate_model_image();

	public:
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

		static void test();

        // Read from an input stream
        void read(std::istream& in);
};

} // namespace Lensing2

#endif

