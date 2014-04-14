#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Lensing2;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* metadata_file)
{
	// Read in the metadata
	fstream fin(metadata_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
	fin>>ni>>nj;
	fin>>x_min>>x_max>>y_min>>y_max;
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
}

