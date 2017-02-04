#ifndef _PSF_
#define _PSF_

#include <armadillo>
#include <ostream>
#include <vector>
#include "DNest4/code/RNG.h"

class PSF
{
	private:
		int size;
		std::vector< std::vector<double> > pixels;

        // Double-gaussian model
        double outer_width;
        double inner_width_frac;
        double inner_mass_frac;
        double q;
        double theta, cos_theta, sin_theta;

		// FFT of the PSF
		arma::cx_mat fft_of_psf;
		bool fft_ready;

	public:
		PSF(int size);

        // DNest4 stuff
        void from_prior(DNest4::RNG& rng);
        double perturb(DNest4::RNG& rng);
        void print(std::ostream& out) const;

        // Get everything ready
        void assemble();

		// Uses loopy method
		void blur_image(std::vector< std::vector<double> >& img) const;
		// Uses fft method
		void blur_image2
			(std::vector< std::vector<double> >& img) const;
		void normalise();
		void calculate_fft(int Ni, int Nj, double psf_power=1.0);
		void load(const char* filename);
		void set_size(int new_size);

		// Getter for pixels
		const std::vector< std::vector<double> >& get_pixels() const
		{ return pixels; }

		// Unit test
		static void test();
};

#endif

