#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Lensing2;

Data Data::instance;

Data::Data()
:psf(1)
,plans_ready(false)
{

}

void Data::load(const char* metadata_file, const char* image_file,
			const char* psf_file)
{
	/*
	* First, read in the metadata
	*/
	int psf_size;
	fstream fin(metadata_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
	fin>>ni>>nj;
	fin>>x_min>>x_max>>y_min>>y_max>>psf_size>>resolution;
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
	image.assign(ni, vector<double>(nj));
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			fin>>image[i][j];
	fin.close();

	/*
	* Load the psf
	*/
	psf.set_size(psf_size);
	psf.load(psf_file);

	/*
	* Initialise the fftw3 plans
	*/
	double* in1 = new double[ni*nj];
	double* in2 = new double[ni*nj];
	fftw_complex* out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(ni/2 + 1)*ni);
	fftw_complex* out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(ni/2 + 1)*ni);
	forwardPlan1 = fftw_plan_dft_r2c_2d(ni, ni,
                                    in1, out1,
                                    FFTW_ESTIMATE);
	forwardPlan2 = fftw_plan_dft_r2c_2d(ni, ni,
				    in2, out2,
				    FFTW_ESTIMATE);
	backPlan = fftw_plan_dft_c2r_2d(ni, ni,
                                    out1, in1,
                                    FFTW_ESTIMATE);
	delete[] in1; delete[] in2;
	fftw_free(out1); fftw_free(out2);
	plans_ready = true;
}

Data::~Data()
{
	if(plans_ready)
	{
		fftw_destroy_plan(forwardPlan1);
		fftw_destroy_plan(forwardPlan2);
		fftw_destroy_plan(backPlan);
	}
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

