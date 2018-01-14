/*
 Double Precision CG - OpenACC
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

double cpu_time_used;
double start_time, end_time;

#include "laplace.cpp"

double ddot(double* v0, double* v1, const int& N){
  double result = 0.0; int i;
  #pragma acc kernels copyin(v0[0:N], v1[0:N])
  {
  #pragma acc loop reduction(+:result)
  for(i=0;i<N;i++) result+=v0[i]*v1[i];
  }
  return result;
}

void daxpyz(double* v0 , double* v1, const double& c , double* v2 , const int& N)
{
  int i;
  #pragma acc kernels copyin(v1[0:N], v2[0:N]), create(v0[0:N])
  {
  for(i=0;i<N;i++) v0[i] = v1[i]+c*v2[i];
  }
  return;
}

void laplacian_times_vector(double* Ax, const double* x, int* site_ip, int* site_im, int* site_jp, int* site_jm, const int& N )
{
  double mass = 0.1; int i;
  #pragma acc data
  {
  #pragma acc kernels loop copyin(x[0:N], site_ip[0:N], site_im[0:N], site_jp[0:N], site_jm[0:N])
  for(i= 0;i<N;i++)
    Ax[i] = (mass*mass-4.0)*x[i] + x[site_ip[i]] + x[site_im[i]] + x[site_jp[i]] + x[site_jm[i]];
  }
  return;
}

void  cg_double(double*, double*, double*, int*, int*, int*, int*, const int&, const double&);

int main(){
  const int L = 1000;          //  lattice size: L x L
  const int N = L*L;
  const double eps_dbl = 0.5;
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

  double* p_dbl = (double*) calloc(N,sizeof(double));
  double* x_dbl = (double*) calloc(N,sizeof(double));
  double* r_dbl = (double*) calloc(N,sizeof(double));

  siteindex(site_ip, site_im, site_jp, site_jm, L);   // set up the indices for neighboring sites on a 2d lattice 

  for( unsigned int i = 0 ; i< N ; i++ ) b[i] = 0.0;
  b[0] = 1.0;                                         // point source at the origin
  printf("matrix size %d \n",L);
  // Set initial vector x_0
  for(unsigned int i=0;i<N;i++) x[i] = 0.0;

  for(unsigned int i=0;i<N;i++) r[i] = b[i];

  bb = ddot(b,b,N);
  rr = ddot(r,r,N);                //  (r,r)
  cr = sqrt(rr/bb);                        //  |r|/|b|

  for(unsigned int i=0;i<N;i++) p_dbl[i] = r[i];

  int iteration = 0;
  start_time = omp_get_wtime();
  while(cr>eps){
    for(unsigned int i=0;i<N;i++) r_dbl[i] = r[i];
    for(unsigned int i=0;i<N;i++) p_dbl[i] = r_dbl[i];   // defect correction

    cg_double(x_dbl, r_dbl, p_dbl, site_ip, site_im, site_jp, site_jm, N, eps_dbl);

    daxpyz(x, x, 1.0, x_dbl, N);              //  x = x + x_dbl

    laplacian_times_vector(Ap, x, site_ip, site_im, site_jp, site_jm, L);

    daxpyz(r, b, -1.0, Ap, N);               //  r = b - A.x

    rr=ddot(r, r, N);                //  (r,r)

    cr = sqrt(rr/bb);                        //  |r|/|b|
    iteration++;
  }
  end_time = omp_get_wtime();

  printf( "=== CG GPU Double Precision has converged ===\n");
  printf( "    Iter,  |r|/|b|: %7d  %16.8E \n", iteration, cr );

  cpu_time_used = ((double) (end_time - start_time));
  printf( "CPU time used %16.8E seconds \n", cpu_time_used);

  return 0;
}

void  cg_double(double* x, double* b, double* p, int* site_ip, int* site_im, int* site_jp, int* site_jm, const int& N, const double& eps)
{

  double* r = new double[N];
  double* Ap = new double[N];        //  A.p
  double* Cp = new double[N];        //  C.p
  double rr;                        //  (r,r)_k
  double tt;                        //  (r,r)_k+1
  double pAp;                       //  (p, A.p)
  double bb;                        //  (b, b)
  double alpha, beta;
  double cr;

// initialization
  for(unsigned int i=0;i<N;i++) x[i] = 0.0;
  for(unsigned int i=0;i<N;i++) r[i] = b[i];

  rr = ddot(r, r, N);
  bb = ddot(b, b, N);
  cr = sqrt(rr/bb);

  unsigned int k= 0;
  while(cr>eps){
    laplacian_times_vector(Ap, p, site_ip, site_im, site_jp, site_jm, N );      //   A.p
    pAp = ddot(p,Ap,N);                                               //   (p, Ap)
    alpha = (rr/pAp);
    daxpyz(x,  x,  alpha, p, N );         //  x = x + alpha * p
    daxpyz(r,  r, -alpha, Ap, N );        //  r = r - alpha * Ap
    tt = ddot(r,r,N);             //  tt = (r, r)
    beta = (tt/rr);
    daxpyz(p, r, beta, p, N);             //  p = r - beta * p
    k++;
    rr = tt;
    cr = sqrt(rr/bb);
  }
  return;
}
