#include <iostream>

using namespace std;

extern "C"
{
	void fastelldefl_(double* x, double* y, double* b, double* gam, double* q,
						double* rcsq, double alpha[]);
}

int main()
{
	double b = 1.; double gam = 1.3;
	double q = 0.9999; double rcsq = 1E-7;
	double alpha[2];
	double x = 1.;
	double y = 0.;

	fastelldefl_(&x, &y, &b, &gam, &q, &rcsq, alpha);
	alpha[0] /= 16.1166;
	alpha[1] /= 16.1166;
	cout<<alpha[0]<<' '<<alpha[1]<<endl;

	return 0;
}

