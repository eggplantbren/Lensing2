#ifndef Lensing2_Source
#define Lensing2_Source

#include <ostream>
#include <vector>
#include "DNest4/code/RNG.h"

namespace Lensing2
{

class Source
{
    protected:


    public:
        virtual ~Source() { }

        // Evaluate at an array of positions
        virtual void evaluate(const std::vector<std::vector<double>>& x,
                              const std::vector<std::vector<double>>& y,
                              std::vector<std::vector<double>>& f,
                              bool update=false) const = 0;

        // MCMC related stuff
        virtual void from_prior(DNest4::RNG& rng) = 0;
        virtual double perturb(DNest4::RNG& rng) = 0;

        virtual void print(std::ostream& out) const = 0;
};

} // namespace Lensing2

#endif

