#ifndef Lensing2_Blobby_h
#define Lensing2_Blobby_h

#include "Source.h"

namespace Lensing2
{

class Blobby:public Source
{
	private:


	public:
		// Required methods
		double evaluate(double x, double y) const;
		void from_prior();
		double perturb();
		void print(std::ostream& out) const;
};

} // namespace Lensing2

#endif

