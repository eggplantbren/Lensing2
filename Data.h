#ifndef _Data_
#define _Data_

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

	public:
		Data();
		void load(const char* metadata_file);

		int get_ni() const { return ni; }
		int get_nj() const { return nj; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

} // namespace Lensing2

#endif

