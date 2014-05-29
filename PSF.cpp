#include "PSF.h"
#include "Utils.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <fftw3.h>

using namespace std;
using namespace DNest3;

PSF::PSF(int size)
:size(size)
,pixels(size, vector<double>(size, 0.))
{
	assert(size%2 == 1);
	pixels[size/2][size/2] = 1.;
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

void PSF::blur_image_using_fftw(vector< vector<double> >& img) const
{
	int Ni = img.size();
	int Nj = img[0].size();
	int ni = pixels.size();
	int nj = pixels[0].size();

	// Make the psf the same size as the image
	vector< vector<double> > psf(Ni, vector<double>(Nj, 0.));
	int m, n;
	for(int i=0; i<ni; i++)
	{
		m = mod(i - ni/2, Ni);
		for(int j=0; j<nj; j++)
		{
			n = mod(j - nj/2, Nj);
			psf[m][n] = pixels[i][j];
		}
	}
}
/*
	// Do the FFT of this one and the other
	fftw_complex* out1; fftw_complex* out2;
	fftw_plan forwardPlan1, forwardPlan2;

	double* in1 = new double[ni*nj];
	double* in2 = new double[ni*nj];

	int k = 0;
	for(unsigned int i=0; i<ni; i++)
		for(unsigned int j=0; j<nj; j++)
		{
			in1[k] = img[i][j];
			in2[k++] = psf[i][j];
		}

	out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(ni/2 + 1)*ni);
	out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(ni/2 + 1)*ni);

	forwardPlan1 = fftw_plan_dft_r2c_2d(ni, ni,
                                    in1, out1,
                                    FFTW_ESTIMATE);
	forwardPlan2 = fftw_plan_dft_r2c_2d(ni, ni,
				    in2, out2,
				    FFTW_ESTIMATE);

	fftw_execute(forwardPlan1);
	fftw_execute(forwardPlan2);

	// (a + bi)(c + di) = ac + adi + bic + bdi^2
	// = (ac - bd) + (ad + bc)i
	// Multiply the ffts and put the result in out1
	double re, im;
	for(unsigned int i=0; i<(ni/2 + 1)*ni; i++)
	{
		re = out1[i][0]*out2[i][0] - out1[i][1]*out2[i][1];
		im = out1[i][0]*out2[i][1] + out2[i][0]*out1[i][1];
		out1[i][0] = re;
		out1[i][1] = im;
	}
	
	delete[] in2;
        fftw_destroy_plan(forwardPlan1);
	fftw_destroy_plan(forwardPlan2);
	fftw_free(out2);

	fftw_plan backPlan = fftw_plan_dft_c2r_2d(ni, ni,
                                    out1, in1,
                                    FFTW_ESTIMATE);

	fftw_execute(backPlan);	

	double coeff = 1.0/(ni*nj);
	k = 0;
	for(unsigned int i=0; i<ni; i++)
		for(unsigned int j=0; j<nj; j++)
			values[nj*i + j] = in1[k++]*coeff;
	fftw_destroy_plan(backPlan);
	delete[] in1;
	fftw_free(out1);
}
*/
