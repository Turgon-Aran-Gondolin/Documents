#include <random>
#include <math.h>
#include "ranlxd.h"


void gauss_rand(int n, double* rand){
	/* std::random_device rd;  */
	/* std::mt19937 gen(rd); */
	double tmp[2];

	std::normal_distribution<double> gau(0,1);

	for (size_t i=0;i<n;i++) {
		/* rand[i] = gau(gen); */
		ranlxd(tmp,2);
		rand[i] = gau(tmp);
	}
}
