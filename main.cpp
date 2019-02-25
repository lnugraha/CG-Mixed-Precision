#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>

#include "libmath.hpp"
#include "laplacian.hpp"

int main(int argc, char *argv[]){

	const unsigned int N = 10;
	single *sin;
	dual *dbl;
	mutual *mut;

	float *vector_0, *vector_1;
	float scalar;
	vector_0 = (float*)calloc(N, sizeof(float));
	vector_1 = (float*)calloc(N, sizeof(float));

	for (unsigned int i=0; i<N; ++i){
		vector_0[i] = 1.5;
		vector_1[i] = 2.0;
	}

	scalar = sin->sdot(vector_0, vector_1, N);
	std::cout << scalar << std::endl;

	free(vector_0);
	free(vector_1);
	return 0;
}
