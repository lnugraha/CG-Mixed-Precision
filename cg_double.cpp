/* Conjugate Gradient (CG) algorithm 
 * To find the solution x of the linear system,  A.x = b
 * A: positive-definite matrix, A = C^\dagger.C
 * b: source vector
 * x: solution vector
 *
 * Initial vector: x_0 (any choice is fine, here x_0 = 0)
 * r_0 = b - A.x_0;
 * p_0 = r_0;
 * bb = (b,b);
 * rr = (r_0, r_0);
 * cr = sqrt(rr/bb); 
 * k = 0; 
 * while (cr > eps) {
 *  Apk = A.p_k;
 *  pAp = (p_k, Apk);
 *  r_(k+1) = r_k - (rr/pAp) Apk;
 *  x_(k+1) = x_k + (rr/pAp) p_k
 *  tt = (r_(k+1), r_(k+1));
 *  p_(k+1) = r_(k+1) + (tt/rr) p_k;
 *  k++;
 *  rr = tt;
 *  cr = sqrt(rr/bb);
 *  } 
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
using namespace std;

clock_t start, end;        // Set up the clock for timing 
double cpu_time_used;
     
#include "Laplace.hpp"
#include "MatrixOperations.hpp"

static const unsigned long long int L = 100;
static const unsigned long long int N = L * L; // 2-D Lattice
static const double eps = 1.0E-06;

int main(int argc, char *argv[]){ 
  unsigned int* site_ip= new unsigned int[N];
  unsigned int* site_im= new unsigned int[N];
  unsigned int* site_jp= new unsigned int[N];
  unsigned int* site_jm= new unsigned int[N];
  
  double* b  = new double[N];
  double* x  = new double[N];
  double* p  = new double[N];
  double* Ap = new double[N];       //  A.p
  double* r  = new double[N];       //  residual vector in double precision
  
  printf("=== To solve A.x = b, using CG with double precision ===\n\n");
  start = clock();   // start the clock
  
  SiteIndices mainSITES;
  mainSITES.site_ip = site_ip;
  mainSITES.site_im = site_im;
  mainSITES.site_jp = site_jp;
  mainSITES.site_jm = site_jm;
  mainSITES.N       = N;


  // siteindex(site_ip, site_im, site_jp, site_jm, L);
  siteindex(mainSITES, L);
  // set up the indices for neighboring sites for 2d lattice

  for (unsigned int i=0 ; i<N ; ++i) b[i] = 0.0;
  b[0] = 1.0; // point source at the origin
  double bb = ddot(b, b, N); // <b, b>
  printf("    b is a source vector of size %d \n\n", N);
  // Set initial vector x_0 
  for (unsigned int i=0; i<N; ++i) { x[i] = 0.0; r[i] = b[i]; p[i] = r[i]; }

  double rr = ddot(r, r, N);  // <r,r>_k
  double cr = sqrt(rr/bb);    // sqrt(rr/bb)
  unsigned int k = 0;
  printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr ); 
  
  dpLaplacianVector DPLP;
  while (cr>eps)
  {
    DPLP.laplacian_times_vector(Ap, p, mainSITES);

    double pAp = ddot(p, Ap, mainSITES.N);  // <p, Ap>
    double alpha = (rr/pAp);
    daxpyz(x, x,  alpha, p, mainSITES.N);   //  x = x + alpha * p
    daxpyz(r, r, -alpha, Ap, mainSITES.N);  //  r = r - alpha * A.p
    double tt = ddot(r, r, mainSITES.N);
    double beta = (tt/rr);
    daxpyz(p, r, beta, p, mainSITES.N);     //  p = r + beta * p
    k++;
    if( (k % 100) == 0 )
    {
    //  printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr );
      printf("    k,  pAp, rr, tt, bb : %10d  %16.8E  %16.8E  %16.8E  %16.8E \n"
      , k, pAp, rr, tt, bb);
    }
    rr = tt;
    cr = sqrt(rr/bb);
  } // END-WHILE

  printf("\n >>>>  CG has converged for eps =%8.1E  <<<< \n",eps);
  printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr ); 

  start = clock() - start;
  cpu_time_used = ((double)(start)) / CLOCKS_PER_SEC;
  printf("\n --- CPU time = %16.8E seconds \n",cpu_time_used);
  
  delete site_ip; delete site_im; delete site_jp; delete site_jm;
  delete b; delete x; delete p; delete r; delete Ap;
  return 0;
}
