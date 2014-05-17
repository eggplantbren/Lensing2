#include "PSF.h"
#include <cassert>

using namespace std;

PSF::PSF(int size)
:pixels(size, vector<double>(size))
{
	assert(size%2 == 1);
}


