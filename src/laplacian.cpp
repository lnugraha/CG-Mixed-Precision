#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "laplacian.hpp"
#include "libmath.hpp"

void laplace::siteindex(int *site_ip, int *site_im,
						int *site_jp, int *site_jm, const int &L){
	unsigned int i,j;
	int  ip, im, jp, jm, site;
	for (i=0; i<L; ++i){
		im = (i-1+L)%L;
		ip = (i+1)%L;
		for (j=0; j<L; ++j){
			site = i*L + j;
			jp   = (j+1)%L;
			jm   = (j-1+L)%L;

			site_ip[site] = ip*L+j;
			site_im[site] = im*L+j;
			site_jp[site] = i*L+jp;
			site_jm[site] = i*L+jm
		} // END-FOR j
	} // END-FOR i

};

void laplace::laplacian_times_vector( double* Ax, const double* x, int* site_ip, int* site_im, 
                             int* site_jp, int* site_jm, const int& N )
{
	double mass = 0.1;
	for(unsigned int i=0; i<N; ++i) Ax[i] = (mass*mass-4.0)*x[i] + x[site_ip[i]] + x[site_im[i]] + x[site_jp[i]] + x[site_jm[i]];
}

void laplace::laplacian_times_vector_single( float* Ax, const float* x, int* site_ip, int* site_im, 
                                    int* site_jp, int* site_jm, const int& N )
{
	float mass = 0.1;
	for(unsigned int i= 0; i<N; ++i) Ax[i] = (mass*mass-4.0)*x[i] + x[site_ip[i]] + x[site_im[i]] + x[site_jp[i]] + x[site_jm[i]];
}
