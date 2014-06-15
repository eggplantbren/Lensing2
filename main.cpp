#include <iostream>
#include <fstream>
#include "RandomNumberGenerator.h"
#include "Start.h"
#include "MyModel.h"
#include "Sources/Blobby.h"
#include "Data.h"

using namespace std;
using namespace DNest3;
using namespace Lensing2;

int main(int argc, char** argv)
{
	// Load some "data"
	Data::get_instance().load("Data/test_metadata.txt",
					"Data/harder_image.txt",
					"Data/test_sigma.txt",
					"Data/test_psf.txt");

	// Run DNest
	MTSampler<MyModel> sampler = setup_mt<MyModel>(argc, argv);
	sampler.run();
	return 0;
}


