#include "Blobby.h"
#include "DNest4/code/RNG.h"

using namespace std;
using namespace DNest4;

namespace Lensing2
{

Blobby::Blobby(double x_min, double x_max,
                    double y_min, double y_max)
:blobs(4, 500, false,
    BasicCircular(x_min, x_max, y_min, y_max), PriorType::log_uniform)
{

}

void Blobby::evaluate(const std::vector<std::vector<double>>& x,
                      const std::vector<std::vector<double>>& y,
                      std::vector<std::vector<double>>& f,
                      bool update) const
{
    // Start the evaluations at zero
    for(size_t i=0; i<f.size(); ++i)
        for(size_t j=0; j<f[i].size(); ++j)
            f[i][j] = 0.0;

    const vector< vector<double> >& components = (update)?
                                                 (blobs.get_added()):
                                                 (blobs.get_components());

    double rsq, widthsq, one_over;
    double c = 2.0/M_PI;

    for(size_t k=0; k<components.size(); ++k)
    {
        widthsq = components[k][3]*components[k][3];
        one_over = 1.0/widthsq;

        for(size_t i=0; i<f.size(); ++i)
        {
            for(size_t j=0; j<f[i].size(); ++j)
            {
                rsq = pow(x[i][j] - components[k][0], 2) +
                                    pow(y[i][j] - components[k][1], 2);

                if(rsq < widthsq)
                    f[i][j] += components[k][2]*c*(1.0 - rsq*one_over)*one_over;
            }
        }
    }
}

void Blobby::from_prior(RNG& rng)
{
    blobs.from_prior(rng);
}

double Blobby::perturb(RNG& rng)
{
    double logH = 0.;

    logH += blobs.perturb(rng);

    return logH;
}

void Blobby::print(ostream& out) const
{
    blobs.print(out);
}

void Blobby::read(istream& in)
{
    blobs.read(in);
}

} // namespace Lensing2

