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
		void normalise();
		void load(const char* filename);
		void set_size(int size);

		// Getter for pixels
		const std::vector< std::vector<double> >& get_pixels() const
		{ return pixels; }

};

#endif
