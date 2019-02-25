#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "libmath.hpp"

float single::sdot(float *v0, float *v1, const int &N){
	float scalar = 0.0;
	for (unsigned int i=0; i<N; ++i) scalar+=v0[i]*v1[i];
	return scalar;
};

void single::saxpyz(float *v0, float *v1, float *v2, const float &c, const int &N){
	for (unsigned int i=0; i<N; ++i) v0[i] = v1[i] + c*v2[i];
};

double dual::ddot(double *v0, double *v1, const int &N){
	double scalar = 0.0;
	for (unsigned int i=0; i<N; ++i) scalar+=v0[i]*v1[i];
	return scalar;
}

void dual::daxpyz(double *v0, double *v1, double *v2, const double &c, const int &N){
	for (unsigned int i=0; i<N; ++i) v0[i] = v1[i] + c*v2[i];
};

void mutual::vec_add(double *vec_dbl, const float *vec_sgl, const int &N){
	for (unsigned int i=0; i<N; ++i) vec_dbl[i] += (double)vec_sgl[i];
};

