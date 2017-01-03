#include <iostream>
#include <fstream>
#include <ctime>
#include "DNest4/code/RNG.h"
#include "Data.h"
#include "BasicCircular.h"

int main()
{
    // Load a dataset
    Lensing2::Data::get_instance().load("../Data/mock_metadata.txt",
                                        "../Data/mock_image.txt",
                                        "../Data/mock_sigma.txt",
                                        "../Data/mock_psf.txt");

    // An RNG
    DNest4::RNG rng(time(0));

    // Something to test
    BasicCircular bc(-1.0, 1.0, -1.0, 1.0);

    for(int i=0; i<1000; ++i)
    {
        bc.from_prior(rng);

        std::vector<double> u(4);
        for(double& uu: u)
            uu = rng.rand();

        std::cout<<"These should be the same:\n";
        for(double& uu: u)
            std::cout<<uu<<' ';
        std::cout<<std::endl;

        bc.from_uniform(u);
        bc.to_uniform(u);
        for(double& uu: u)
            std::cout<<uu<<' ';
        std::cout<<'\n'<<std::endl;
    }

    // Now use the Metropolis algorithm
    bc.from_prior(rng);
    bc.print(std::cout); std::cout<<'\n';
    std::vector<double> x{0.5, 0.5, 0.5, 0.5};
    bc.from_uniform(x);
    double logp = bc.log_pdf(x);
    std::fstream fout("output.txt", std::ios::out);
    std::cout<<"Now doing metropolis..."<<std::endl;

    for(int i=0; i<10000000; ++i)
    {
        auto proposal = x;
        int which = rng.rand_int(proposal.size());

        proposal[which] += exp(5*rng.randn())*rng.randn();
        double logp2 = bc.log_pdf(proposal);

        if(rng.rand() <= exp(logp2 - logp))
        {
            x = proposal;
            logp = logp2;
        }

        if((i+1) % 10 == 0)
        {
            for(double xx: x)
                fout<<xx<<' ';
            fout<<std::endl;
        }
    }
    fout.close();

    return 0;
}

