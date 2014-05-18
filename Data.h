#ifndef _Data_
#define _Data_

#include <vector>
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

		// The ray grid
		std::vector< std::vector<double> > x_rays;
		std::vector< std::vector<double> > y_rays;

		// The pixels
		std::vector< std::vector<double> > image;

		// The PSF
		PSF psf;

		void compute_ray_grid();

	public:
		Data();
		void load(const char* metadata_file, const char* image_file,
				const char* psf_file);

		// Getters
		int get_ni() const { return ni; }
		int get_nj() const { return nj; }
		double get_x_min() const { return x_min; }
		double get_x_max() const { return x_max; }
		double get_y_min() const { return y_min; }
		double get_y_max() const { return y_max; }
		const std::vector< std::vector<double> >& get_x_rays() const
			{ return x_rays; }
		const std::vector< std::vector<double> >& get_y_rays() const
			{ return y_rays; }
		const std::vector< std::vector<double> >& get_image() const
			{ return image; }
		const PSF& get_psf() const { return psf; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

} // namespace Lensing2

#endif

