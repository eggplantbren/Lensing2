#include <iostream>
#include <fstream>
#include "RandomNumberGenerator.h"
#include "Sources/Blobby.h"
#include "Data.h"

using namespace std;
using namespace DNest3;
using namespace Lensing2;

void test1();
void test2();

int main()
{
	test1();
	return 0;
}

void test1()
{
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(2);

	Blobby b(-1., 1., -1., 1., 1E-3, 1E3);
	b.from_prior();
	Blobby b2 = b;

	int steps = 10000;
	int skip = 10;

	fstream fout("output.txt", ios::out);
	// Explore the prior for Blobbys
	for(int i=0; i<steps; i++)
	{
		b2 = b;
		double logH = b2.perturb();
		if(randomU() <= exp(logH))
			b = b2;

		if(i%skip == 0)
		{
			for(double y = 0.99; y >= -0.99; y -= 0.01)
			{
				for(double x = -0.99; x <= 0.99; x += 0.01)
					fout<<b.evaluate(x, y)<<' ';
			}
			fout<<endl;
			cout<<i<<endl;
		}
	}

	fout.close();
}

void test2()
{
	Data::get_instance().load("Data/test_metadata.txt");
}

