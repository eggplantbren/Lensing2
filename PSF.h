#ifndef _PSF_
#define _PSF_

#include <vector>

class PSF
{
	private:
		int size;
		std::vector< std::vector<double> > pixels;

	public:
		PSF(int size);
		void blur_image(std::vector< std::vector<double> >& img) const;
//		void blur_image_using_fftw
//			(std::vector< std::vector<double> >& img) const;
		void normalise();
		void load(const char* filename);
		void set_size(int new_size);

		// Getter for pixels
		const std::vector< std::vector<double> >& get_pixels() const
		{ return pixels; }

		// Unit test
		static void test();

};


#endif

