#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Lensing2;

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
	int psf_size;
	fstream fin(metadata_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
	fin>>ni>>nj;
	fin>>x_min>>x_max>>y_min>>y_max>>psf_size>>resolution>>fft_flag2;
	fin.close();

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
	fin.open(image_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<image_file<<"."<<endl;
	image = arma::mat(ni, nj);
	for(int i=0; i<ni; ++i)
		for(int j=0; j<nj; ++j)
			fin>>image(i, j);
	fin.close();

	/*
	* Load the sigma map
	*/
	fin.open(sigma_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<sigma_file<<"."<<endl;
	sigma = arma::mat(ni, nj);
	for(int i=0; i<ni; ++i)
		for(int j=0; j<nj; ++j)
			fin>>sigma(i, j);
	fin.close();

	/*
	* Load the psf
	*/
	psf.set_size(psf_size);
	psf.load(psf_file);

//	if(psf_size >= 15)
//	{
		fft_flag1 = true;
		if(fft_flag2)
			psf.calculate_fft(ni*resolution, nj*resolution);
		else
			psf.calculate_fft(ni, nj);
//	}
}

void Data::compute_ray_grid()
{
	// Make vectors of the correct size
	x_rays = arma::mat(ni*resolution, nj*resolution);
	y_rays = arma::mat(ni*resolution, nj*resolution);

	// Distance between adjacent rays
	double L = dx/resolution;

    for(size_t j=0; j<x_rays.n_cols; ++j)
    {
        for(size_t i=0; i<x_rays.n_rows; ++i)
		{
			x_rays(i, j) = x_min + (j + 0.5)*L;
			y_rays(i, j) = y_max - (i + 0.5)*L;
		}
	}
}

