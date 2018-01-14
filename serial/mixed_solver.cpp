/*
 Mixed Precision CG

   while ( |r|/|b| > eps ) {
     r := b - A.x
     bs := r
     Solve As.xs = bs in single precision
     x := x + xs
   }
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

clock_t start,end;
double cpu_time_used;

#include "laplace.cpp"
#include "axpy.cpp"
#include "matrix_vector.cpp"
#include "cg_single.cpp"
#include "vector_addition_dbl_sgl.cpp"

int main(){
  const int L = 8000;          //  lattice size: L x L
  const int N = L*L;
  const float eps_sgl = 0.25;
  const double eps = 1.0e-6;

  int* site_ip = (int*) calloc(N,sizeof(int));
  int* site_im = (int*) calloc(N,sizeof(int));
  int* site_jp = (int*) calloc(N,sizeof(int));
  int* site_jm = (int*) calloc(N,sizeof(int));

  double* b  = (double*) calloc(N,sizeof(double));
  double* x  = (double*) calloc(N,sizeof(double));
  double* p  = (double*) calloc(N,sizeof(double));
  double* Ap = (double*) calloc(N,sizeof(double));       //  A.p
  double* r  = (double*) calloc(N,sizeof(double));;        //  residual vector in double precision
  double bb;                        //  (b,b)
  double rr;                        //  (r,r)
  double cr;                        //  sqrt(rr/bb)
  float* r_sgl = (float*) calloc(N,sizeof(float));      // residual vector in single precision
  float* x_sgl = (float*) calloc(N,sizeof(float));
  float* p_sgl = (float*) calloc(N,sizeof(float));
  float  rs;                        //  sqrt(r_sgl,r_sgl)

  siteindex(site_ip, site_im, site_jp, site_jm, L);   // set up the indices for neighboring sites on a 2d lattice 

  for(unsigned int i = 0 ; i< N ; i++ ) b[i] = 0.0;
  b[0] = 1.0;                                         // point source at the origin
  printf("matrix size %d \n",L);

  for(unsigned int i=0;i<N;i++) x[i] = 0.0;
  for(unsigned int i=0;i<N;i++) r[i] = b[i];

  bb = ddot(b, b, N);
  rr = ddot(r, r, N);                //  (r,r)
  cr = sqrt(rr/bb);                        //  |r|/|b|

  for(unsigned int i=0;i<N;i++) p_sgl[i] = (float) r[i];

  int iteration = 0;
  start = clock();
  while(cr>eps) {
    for(unsigned int i=0;i<N;i++) r_sgl[i] = (float) r[i];
    for(unsigned int i=0;i<N;i++) p_sgl[i] = r_sgl[i];   // defect correction

    cg_single(x_sgl, r_sgl, p_sgl, site_ip, site_im, site_jp, site_jm, N, eps_sgl);
    vector_addition_dbl_sgl(x, x_sgl, N);

    laplacian_times_vector(Ap, x, site_ip, site_im, site_jp, site_jm, L);
    daxpyz(r, b, -1.0, Ap, N);               //  r = b - A.x
    rr = ddot(r, r, N);                //  (r,r)
    cr = sqrt(rr/bb);                        //  |r|/|b|
    iteration++;
  }
  end = clock();
  printf( "=== CG with Mixed Precision has converged ===\n");
  printf( "    Iter,  |r|/|b|: %7d  %16.8E \n", iteration, cr );

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf( "CPU time used %16.8E seconds \n",cpu_time_used);

  return 0;
}
