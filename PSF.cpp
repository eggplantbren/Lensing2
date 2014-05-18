#include "PSF.h"
#include <cassert>

using namespace std;

PSF::PSF(int size)
:pixels(size, vector<double>(size))
{
	assert(size%2 == 1);
}

void PSF::blur_image(vector< vector<double> >& img) const
{
	// Make result image. Assume img is rectangular...
	vector< vector<double> > result(img.size(),
					vector<double>(img[0].size(), 0.));
}

