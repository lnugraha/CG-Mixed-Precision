#include "MatrixOperations.hpp"

double ddot( const double* w, const double* v, const unsigned int& N )
{
  double norm = 0.0;
  for (unsigned int i=0; i<N; ++i) norm += w[i] * v[i];
  return norm;
}

float sdot( const float* w, const float* v, const unsigned int& N )
{
  float norm = 0.0;
  for (unsigned int i=0; i<N; ++i) norm += w[i] * v[i];
  return norm;
}

void saxpyz( float* v0, float* v1, const float& c, float* v2,
  const unsigned int& N )
{ for (unsigned int i=0; i<N; i++) v0[i] = v1[i] + c * v2[i]; }

void daxpyz( double* v0, double* v1, const double& c, double* v2, 
  const unsigned int& N )
{ for (unsigned int i=0; i<N; i++) v0[i] = v1[i] + c * v2[i]; }

void vector_addition_dbl_sgl( double* x_dbl, const float* x_sgl, 
  const unsigned int& N )
{ for (unsigned int i=0; i<N; i++) x_dbl[i] += (double)x_sgl[i]; }

// ========================================================================== //

void dpLaplacianVector::laplacian_times_vector( double* Ax, const double* x, 
  const SiteIndices& SITES, const unsigned int& limit)
{
  static double mass = 0.001;
  for (unsigned int i=0; i<limit; i++)
  {
    Ax[i] = (mass*mass-4.0) * x[i] + x[ SITES.site_ip[i] ] + 
    x[ SITES.site_im[i] ] + x[ SITES.site_jp[i] ] + x[ SITES.site_jm[i] ];
  }
}

void spLaplacianVector::laplacian_times_vector( float* Ax, const float* x, 
  const SiteIndices& SITES, const unsigned int& limit)
{
  static float mass = 0.001;
  for (unsigned int i=0; i<limit; i++)
  {
    Ax[i] = (mass*mass-4.0) * x[i] + x[ SITES.site_ip[i] ] + 
    x[ SITES.site_im[i] ] + x[ SITES.site_jp[i] ] + x[ SITES.site_jm[i] ];
  }
}
