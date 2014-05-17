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
	// Make backup
	vector< vector<double> > backup = img;
}

