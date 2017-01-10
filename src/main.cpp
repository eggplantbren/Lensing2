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
	Data::get_instance().load("Data/horseshoe_metadata.txt",
					"Data/horseshoe_image.txt",
					"Data/horseshoe_sigma.txt",
					"Data/horseshoe_psf.txt");

	// Run DNest
	Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
	sampler.run();
	return 0;
}


