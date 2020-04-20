#ifndef CGSOLVERS_HPP_
#define CGSOLVERS_HPP_

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include "MatrixOperations.hpp"

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
