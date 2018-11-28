#include <iostream>
#include <fstream>
#include "yaml-cpp/yaml.h"
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Sources/Blobby.h"
#include "Data.h"

using namespace std;
using namespace DNest4;
using namespace Lensing2;

int main(int argc, char** argv)
{

    // Load the run information from run.yaml
    YAML::Node config;
    try
    {
        config = YAML::LoadFile("run_files.yaml");
    }
    catch(...)
    {
        std::cerr << "# Couldn't open or parse run_files.yaml." << std::endl;
        std::cerr << "# Aborting." << std::endl;
        return 0;
    }

//    // Read in the values
//    prior_max_power = config["prior"]["max_power"].as<double>();
//    prior_scale = config["prior"]["scale"].as<double>();
//    grid_num_points = config["grid"]["num_points"].as<size_t>();
//    lms_max_power = config["lagrange_multipliers"]["max_power"].as<double>();
//    lms_scale = config["lagrange_multipliers"]["scale"].as<double>();

//    std::cout << "done." << std::endl;
//}

	// Load the data
	Data::get_instance().load(config["metadata_file"].as<std::string>().c_str(),
					          config["image_file"].as<std::string>().c_str(),
					          config["sigma_file"].as<std::string>().c_str(),
					          config["psf_file"].as<std::string>().c_str());

	// Run DNest
    RNG::randh_is_randh2 = true;
	Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
	sampler.run();
	return 0;
}


