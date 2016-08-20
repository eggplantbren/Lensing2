#ifndef _Data_
#define _Data_

#include <vector>
#include <armadillo>
#include "PSF.h"

namespace Lensing2
{

class Data
{
	private:
		// Number of pixels
		int ni, nj;

		// Coordinates of image boundaries
		double x_min, x_max, y_min, y_max;

		// Pixel widths
		double dx, dy;

		// Number of rays per pixel per dimension
		int resolution;

		// The ray grid
		arma::mat x_rays;
		arma::mat y_rays;

		// The pixels
		arma::mat image;

		// Sigma map
		arma::mat sigma;

		// The PSF
		PSF psf;

		// Flag1: use FFTs to do blurring?
		// Flag2: Apply blurring at the higher resolution?
		bool fft_flag1, fft_flag2;

		void compute_ray_grid();

	public:
		Data();
		void load(const char* metadata_file, const char* image_file,
				const char* sigma_file, const char* psf_file);

		// Getters
		int get_ni() const { return ni; }
		int get_nj() const { return nj; }
		double get_x_min() const { return x_min; }
		double get_x_max() const { return x_max; }
		double get_y_min() const { return y_min; }
		double get_y_max() const { return y_max; }
		int get_resolution() const { return resolution; }
		const arma::mat& get_x_rays() const
			{ return x_rays; }
		const arma::mat& get_y_rays() const
			{ return y_rays; }
		const arma::mat& get_image() const
			{ return image; }
		const arma::mat& get_sigma() const
			{ return sigma; }
		const PSF& get_psf() const { return psf; }
		bool use_fft() const { return fft_flag1; }
		bool psf_is_highres() const { return fft_flag2; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

} // namespace Lensing2

#endif

