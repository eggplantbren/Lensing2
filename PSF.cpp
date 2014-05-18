#include "PSF.h"
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

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
			if(img[i][j] > 0.)
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

