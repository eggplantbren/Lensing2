#include <iostream>
#include <fstream>
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Sources/Blobby.h"
#include "Data.h"

using namespace std;
using namespace DNest4;
using namespace Lensing2;

int main(int argc, char** argv)
{
	// Load some "data"
	Data::get_instance().load("Data/mock_metadata.txt",
					"Data/harder_image.txt",
					"Data/mock_sigma.txt",
					"Data/mock_psf.txt");

	// Run DNest
	Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
	sampler.run();
	return 0;
}

