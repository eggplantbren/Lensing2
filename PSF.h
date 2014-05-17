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


};

#endif

