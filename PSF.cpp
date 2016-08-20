#include "PSF.h"
#include "DNest4/code/Utils.h"
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;
using namespace DNest4;
using namespace arma;

PSF::PSF(int size)
:size(size)
,pixels(size, size)
,fft_ready(false)
{
	assert(size%2 == 1);
	pixels(size/2, size/2) = 1.;
}

void PSF::set_size(int new_size)
{
	size = new_size;
	pixels = arma::mat(size, size);
	pixels(size/2, size/2) = 1.;
}

void PSF::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# ERROR: couldn't open file "<<filename<<"."<<endl;
		return;
	}
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			fin>>pixels(i, j);
	fin.close();
	normalise();
}

void PSF::normalise()
{
	double sum = 0.;
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			sum += pixels(i, j);

	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			pixels(i, j) /= sum;
}

void PSF::calculate_fft(int Ni, int Nj)
{
	// Make the psf the same size as the image
	mat psf(Ni, Nj);
	psf.zeros();

	int ni = pixels.n_rows;
	int nj = pixels.n_cols;

	int m, n;
	for(int i=0; i<ni; i++)
	{
		m = mod(i - ni/2, Ni);
		for(int j=0; j<nj; j++)
		{
			n = mod(j - nj/2, Nj);
			psf(m, n) = pixels(i, j);
		}
	}

	fft_of_psf = fft2(psf);
	fft_ready = true;
}


void PSF::blur_image2(mat& img) const
{
	if(!fft_ready)
		cerr<<"# Blurring failed."<<endl;

    // A = img
	// Do the fft of it
	cx_mat B = fft2(img);

	// Multiply the two ffts
	for(size_t j=0; j<img.n_cols; ++j)
        for(size_t i=0; i<img.n_rows; ++i)
			B(i, j) *= fft_of_psf(i, j);

	// Do the inverse fft
	B = ifft2(B);

	// Put back in img
	for(size_t j=0; j<img.n_cols; ++j)
        for(size_t i=0; i<img.n_rows; ++i)
			img(i, j) = real(B(i, j));
}


