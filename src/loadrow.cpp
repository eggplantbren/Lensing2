#include <iostream>
#include <fstream>
#include "MyModel.h"
#include "Data.h"

/*
* Loads a model from posterior_sample.txt and prints it to screen again.
*/

int main(int argc, char** argv)
{
    if(argc <= 1)
    {
        std::cerr<<"Usage: loadrow <row_number>"<<std::endl;
        return 0;
    }

	// Load some "data" --- this is needed!
    Lensing2::Data::get_instance().load("Data/mock_metadata.txt",
                            "Data/harder_image.txt",
                            "Data/mock_sigma.txt",
                            "Data/mock_psf.txt");

    // Model object to read into
    Lensing2::MyModel m;

    std::fstream fin("posterior_sample.txt", std::ios::in);
    // Read past comment lines at the top of the file
    while(fin.peek() == '#')
        fin.ignore(10000000, '\n');

    // Skip through to requested row
    for(int i=0; i<atoi(argv[1]); ++i)
        fin.ignore(10000000, '\n');

    m.read(fin);
    fin.close();

    m.print(std::cout); std::cout<<std::endl;

    return 0;
}

