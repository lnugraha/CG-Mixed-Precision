/*
  Conjugate Gradient (CG) algorithm
*/

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
