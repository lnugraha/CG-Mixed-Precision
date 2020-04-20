#include "CGSolvers.hpp"

void dpCG::Solver( double* x, const double* b, double *p,
  const SiteIndices& SITES, const double& eps)
{
  printf( "=== Double Precision ===\n" );
  double* r = new double[SITES.N];
  double* Ap = new double[SITES.N];        //  A.p
  double* Cp = new double[SITES.N];        //  C.p

  // initialization  
  for (unsigned int i=0; i<SITES.N; i++)
  { x[i] = 0.0; r[i] = b[i]; }
  double rr = ddot(r, r, SITES.N);
  double bb = ddot(b, b, SITES.N);
  double cr = sqrt(rr/bb);
  // end of initialization

  unsigned int k= 0; dpLaplacianVector DPLP;
  while (cr>eps)
  {
    DPLP.laplacian_times_vector(Ap, p, SITES, SITES.N);
    double pAp = ddot(p, Ap, SITES.N);
    double alpha = (rr/pAp);
    daxpyz(x,  x,  alpha, p, SITES.N );         //  x = x + alpha * p
    daxpyz(r,  r, -alpha, Ap, SITES.N );        //  r = r - alpha * Ap
    double tt = ddot(r, r, SITES.N);
    double beta = (tt/rr);
    daxpyz(p, r, beta, p, SITES.N);             //  p = r - beta * p
    k++;
    rr = tt;
    cr = sqrt(rr/bb);
    printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr ); 
  }
  delete r; delete Ap; delete Cp;
}

void spCG::Solver( float* x, const float* b, float *p,
  const SiteIndices& SITES, const float& eps)
{
  printf( "=== Single Precision ===\n" );
  float* r  = new float[SITES.N];
  float* Ap = new float[SITES.N];        //  A.p
  float* Cp = new float[SITES.N];        //  C.p

  // initialization  
  for( unsigned int i=0; i<SITES.N; i++)
  { x[i] = 0.0; r[i] = b[i]; }
  float rr = sdot(r, r, SITES.N);
  float bb = sdot(b ,b, SITES.N);
  float cr = sqrt(rr/bb);          
  // end of initialization

  unsigned int k = 0; spLaplacianVector SPLP;
  while (cr>eps)
  {
    SPLP.laplacian_times_vector(Ap, p, SITES, SITES.N);
    float pAp = sdot(p, Ap, SITES.N);
    float alpha = (rr/pAp);
    saxpyz(x, x,  alpha, p, SITES.N );        //  x = x + alpha * p
    saxpyz(r, r, -alpha, Ap, SITES.N );       //  r = r - alpha * Ap
    float tt = sdot(r, r, SITES.N);
    float beta = (tt/rr);
    saxpyz(p, r, beta, p, SITES.N);           //  p = r - beta * p
    k++;
    rr = tt;
    cr = sqrt(rr/bb);
    printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr ); 
  }
  delete r; delete Ap; delete Cp;
}
