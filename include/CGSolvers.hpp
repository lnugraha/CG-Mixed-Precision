#ifndef CGSOLVERS_HPP_
#define CGSOLVERS_HPP_

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include "MatrixOperations.hpp"

void  cg_single( float* x, const float* b, float* p, 
  const unsigned int* site_ip, const unsigned int* site_im, 
  const unsigned int* site_jp, const unsigned int* site_jm, 
  const unsigned int& N, const float& eps);

void  cg_double( double* x, const double* b, double* p, 
  const unsigned int* site_ip, const unsigned int* site_im, 
  const unsigned int* site_jp, const unsigned int* site_jm, 
  const unsigned int& N, const double& eps);

class CG{
  protected:
    SiteIndices SITES;
  public:

};

class dpCG: public CG{
  public:

    void Solver( double* x, const double* b, double *p,
    const SiteIndices& SITES, const double& eps);
};

class spCG: public CG{
  public:

    void Solver( float* x, const float* b, float *p,
    const SiteIndices& SITES, const float& eps);
};

#endif
