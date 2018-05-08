#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "yaml-cpp/yaml.h"

using namespace std;

namespace Lensing2
{

Data Data::instance;

Data::Data()
:psf(1)
,fft_flag1(false)
,fft_flag2(false)
{

}

void Data::load(const char* metadata_file, const char* image_file,
			    const char* sigma_file, const char* psf_file)
{
	/*
	* First, read in the metadata
	*/
    YAML::Node yaml;
    try
    {
        yaml = YAML::LoadFile(metadata_file);
    }
    catch(...)
    {
        std::cerr << "# Couldn't open or parse " << metadata_file << ".";
        std::cerr << std::endl;
        std::cerr << "# Aborting." << std::endl;
        return;
    }

    // Extract the values
    ni = yaml["dimensions"]["ni"].as<int>();
    nj = yaml["dimensions"]["nj"].as<int>();
    x_min = yaml["dimensions"]["x_min"].as<double>();
    x_max = yaml["dimensions"]["x_max"].as<double>();
    y_min = yaml["dimensions"]["y_min"].as<double>();
    y_max = yaml["dimensions"]["y_max"].as<double>();

    int psf_size = yaml["psf"]["num_pixels"].as<int>();
    resolution = yaml["computation"]["nrays"].as<int>();
    fft_flag2 = yaml["psf"]["is_highres"].as<bool>();

	// Make sure maximum > minimum
	if(x_max <= x_min || y_max <= y_min)
		cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<endl;

	// Compute pixel widths
	dx = (x_max - x_min)/nj;
	dy = (y_max - y_min)/ni;

	// Check that pixels are square
	if(abs(log(dx/dy)) >= 1E-3)
		cerr<<"# ERROR: pixels aren't square."<<endl;

	compute_ray_grid();

	/*
	* Now, load the image
	*/
    std::fstream fin;
	fin.open(image_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<image_file<<"."<<endl;
	image.assign(ni, vector<double>(nj));
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			fin>>image[i][j];
	fin.close();

	/*
	* Load the sigma map
	*/
	fin.open(sigma_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<sigma_file<<"."<<endl;
	sigma.assign(ni, vector<double>(nj));
	for(size_t i=0; i<sigma.size(); i++)
		for(size_t j=0; j<sigma[i].size(); j++)
			fin>>sigma[i][j];
	fin.close();

	/*
	* Load the psf
	*/
	psf.set_size(psf_size);
	psf.load(psf_file);

	if(psf_size >= 15)
	{
		fft_flag1 = true;
		if(fft_flag2)
			psf.calculate_fft(ni*resolution, nj*resolution);
		else
			psf.calculate_fft(ni, nj);
	}

    // Now save the filenames used to run_data.txt
    fstream fout("run_data.txt", ios::out);
    fout<<metadata_file<<endl;
    fout<<image_file<<endl;
    fout<<sigma_file<<endl;
    fout<<psf_file<<endl;
    fout.close();
}

void Data::compute_ray_grid()
{
	// Make vectors of the correct size
	x_rays.assign(ni*resolution,
			vector<double>(nj*resolution));
	y_rays.assign(ni*resolution,
			vector<double>(nj*resolution));

	// Distance between adjacent rays
	double L = dx/resolution;

	for(size_t i=0; i<x_rays.size(); i++)
	{
		for(size_t j=0; j<x_rays[i].size(); j++)
		{
			x_rays[i][j] = x_min + (j + 0.5)*L;
			y_rays[i][j] = y_max - (i + 0.5)*L;
		}
	}
}

} // namespace Lensing2

