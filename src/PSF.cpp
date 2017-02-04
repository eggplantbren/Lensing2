#include "PSF.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include "DNest4/code/Utils.h"

using namespace std;
using namespace DNest4;
using namespace arma;

PSF::PSF(int size)
:size(size)
,pixels(size, vector<double>(size, 0.))
,fft_ready(false)
{
	assert(size%2 == 1);
	pixels[size/2][size/2] = 1.;
}

void PSF::from_prior(DNest4::RNG& rng)
{
    outer_width = size*rng.rand();
    inner_width_frac = rng.rand();
    inner_mass_frac = rng.rand();
    q = exp(0.2*rng.randn());
    theta = 2*M_PI*rng.rand();
    cos_theta = cos(theta); sin_theta = sin(theta);
}

double PSF::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(5);
    if(which == 0)
    {
        outer_width += size*rng.randh();
        DNest4::wrap(outer_width, 0.0, size);
    }
    else if(which == 1)
    {
        inner_width_frac += rng.randh();
        DNest4::wrap(inner_width_frac, 0.0, 1.0);
    }
    else if(which == 2)
    {
        inner_mass_frac += rng.randh();
        DNest4::wrap(inner_mass_frac, 0.0, 1.0);
    }
    else if(which == 3)
    {
        q = log(q);
        logH -= -0.5*pow(q/0.2, 2);
        q += 0.2*rng.randh();
        logH += -0.5*pow(q/0.2, 2);
        q = exp(q);
    }
    else
    {
        theta += 2.0*M_PI*rng.randh();
        DNest4::wrap(theta, 0.0, 2*M_PI);
    }

    return logH;
}

void PSF::set_size(int new_size)
{
	size = new_size;
	pixels.assign(size, vector<double>(size, 0.));
	pixels[size/2][size/2] = 1.;
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
			fin>>pixels[i][j];
	fin.close();
	normalise();
}

void PSF::normalise()
{
	double sum = 0.;
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			sum += pixels[i][j];

	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			pixels[i][j] /= sum;
}

void PSF::calculate_fft(int Ni, int Nj, double psf_power)
{
	// Make the psf the same size as the image
	mat psf(Ni, Nj);
	psf.zeros();

	int ni = pixels.size();
	int nj = pixels[0].size();

	int m, n, sign;
	for(int i=0; i<ni; i++)
	{
		m = mod(i - ni/2, Ni);
		for(int j=0; j<nj; j++)
		{
			n = mod(j - nj/2, Nj);

            sign = 1;
            if(pixels[i][j] < 0.0)
                sign = -1;

            psf(m, n) = sign*pow(std::abs(pixels[i][j]), psf_power);
		}
	}

    // Normalise back to 1
    double tot = 0.0;
    for(size_t j=0; j<psf.n_cols; ++j)
        for(size_t i=0; i<psf.n_rows; ++i)
            tot += psf(i, j);
    for(size_t j=0; j<psf.n_cols; ++j)
        for(size_t i=0; i<psf.n_rows; ++i)
        psf(i, j) /= tot;

    fft_of_psf = fft2(psf);
    fft_ready = true;
}

void PSF::blur_image(vector< vector<double> >& img) const
{
	// Make result image. Assume img is rectangular...
	vector< vector<double> > result(img.size(),
					vector<double>(img[0].size(), 0.));

	int h = size/2;
	int ii, jj;
	int M = static_cast<int>(img.size());
	int N = static_cast<int>(img[0].size());

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			if(img[i][j] != 0.)
			{
				for(int m=0; m<size; m++)
				{
					ii = i + m - h;
					for(int n=0; n<size; n++)
					{
						jj = j + n - h;
						if(ii >= 0 && ii < M &&
							jj >= 0 && jj < N)
							result[ii][jj] +=
							img[i][j]*pixels[m][n];
					}
				}
			}
		}
	}

	img = result;
}

void PSF::blur_image2(vector< vector<double> >& img) const
{
	if(!fft_ready)
		cerr<<"# Blurring failed."<<endl;

	// Copy the image into an Armadillo matrix
	mat A(img.size(), img[0].size());
	for(size_t i=0; i<img.size(); i++)	
		for(size_t j=0; j<img[0].size(); j++)
			A(i, j) = img[i][j];

	// Do the fft of it
	cx_mat B = fft2(A);

	// Multiply the two ffts
	for(size_t i=0; i<img.size(); i++)
		for(size_t j=0; j<img[0].size(); j++)
			B(i, j) *= fft_of_psf(i, j);

	// Do the inverse fft
	B = ifft2(B);

	// Put back in img
	for(size_t i=0; i<img.size(); i++)
		for(size_t j=0; j<img[0].size(); j++)
			img[i][j] = real(B(i, j));

}


void PSF::test()
{
	PSF psf(5);
	psf.load("Data/test_psf.txt");
	psf.calculate_fft(20, 20);

	// Point image
	vector< vector<double> > pixels(20, vector<double>(20, 0.));
	pixels[10][10] = 1.;

	psf.blur_image(pixels);
	// Print image and reset it to zero
	for(size_t i=0; i<pixels.size(); i++)
	{
		for(size_t j=0; j<pixels.size(); j++)
		{
			cout<<pixels[i][j]<<' ';
			pixels[i][j] = 0.;
		}
	}
	cout<<endl;

	// Do it again with ffts
	pixels[10][10] = 1.;

	psf.blur_image2(pixels);
	for(size_t i=0; i<pixels.size(); i++)
		for(size_t j=0; j<pixels.size(); j++)
			cout<<pixels[i][j]<<' ';
	cout<<endl;
}


