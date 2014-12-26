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
	Data::get_instance().load("Data/mock_metadata.txt",
					"Data/mock_image.txt",
					"Data/mock_sigma.txt",
					"Data/mock_psf.txt");

	// Run DNest
	MTSampler<MyModel> sampler = setup_mt<MyModel>(argc, argv);
	sampler.run();
	return 0;
}


