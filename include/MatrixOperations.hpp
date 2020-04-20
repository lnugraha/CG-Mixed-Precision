#ifndef MATRIXOPERATION_HPP_
#define MATRIXOPERATION_HPP_

#include <cstdio>
#include <cstdlib>
#include <cmath>

typedef struct SiteIndices{
  unsigned int *site_ip, *site_im, *site_jp, *site_jm;
  unsigned int N;
} SiteIndices;

// TODO: Merge with double precision library
float sdot( const float* w, const float* v, const unsigned int& N);

double ddot( const double* w, const double* v, const unsigned int& N);

void saxpyz( float* v0, float* v1, const float& c, float* v2, 
  const unsigned int& N );

void daxpyz( double* v0, double* v1, const double& c, double* v2, 
  const unsigned int& N );

void vector_addition_dbl_sgl( double* x_dbl, const float* x_sgl, 
  const unsigned int& N );



void laplacian_times_vector( double* Ax, const double* x, 
  const unsigned int* site_ip, const unsigned int* site_im, 
  const unsigned int* site_jp, const unsigned int* site_jm, 
  const unsigned int& N );



void laplacian_times_vector_single( float* Ax, const float* x, 
  const unsigned int* site_ip, const unsigned int* site_im, 
  const unsigned int* site_jp, const unsigned int* site_jm, 
  const unsigned int& N );

class Laplacian_Times_Vector{
  protected:
    SiteIndices SITES;
  public:

};

class dpLaplacianVector : public Laplacian_Times_Vector{
  public:
    void laplacian_times_vector( double* Ax, const double* x, 
    const SiteIndices& SITES );
};

class spLaplacianVector : public Laplacian_Times_Vector{
  public:
    void laplacian_times_vector( float* Ax, const float* x, 
    const SiteIndices& SITES );
};

#endif
