#ifndef Lensing2_Blobby
#define Lensing2_Blobby

#include "../Source.h"
#include "../BasicCircular.h"
#include "DNest4/code/RJObject/RJObject.h"

namespace Lensing2
{

class Blobby:public Source
{
	private:
		DNest4::RJObject<BasicCircular> blobs;

	public:
		// Pass in image scale and flux scale
		Blobby(double x_min, double x_max,
					double y_min, double y_max);

		// Required methods
		double evaluate(double x, double y, bool update) const;
		void from_prior(DNest4::RNG& rng);
		double perturb(DNest4::RNG& rng);
		void print(std::ostream& out) const;

        // Getter
		const DNest4::RJObject<BasicCircular>& get_blobs() const
        { return blobs; }
};

} // namespace Lensing2

#endif

