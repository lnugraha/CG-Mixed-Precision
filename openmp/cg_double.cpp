/*
  Conjugate Gradient (CG) algorithm

  To find the solution x of the linear system,  A.x = b

  A: positive-definite matrix, A = C^\dagger.C
  b: source vector
  x: solution vector

 Initial vector: x_0 (any choice is fine, here we set x_0 = 0)

  r_0 = b - A.x_0;
  p_0 = r_0;
  bb = (b,b);
  rr = ( r_0, r_0 );
  cr = sqrt(rr/bb);
  k = 0;
  while ( cr > eps ) {
    Apk = A.p_k;
    pAp = (p_k, Apk);
    x_(k+1) = x_k + (rr/pAp) p_k
    r_(k+1) = r_k - (rr/pAp) A.p_k;
    tt = (r_(k+1), r_(k+1));
    p_(k+1) = r_(k+1) + (tt/rr) p_k;
    k++;
    rr = tt;
    cr = sqrt(rr/bb);
  }

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

// initialization
  for(unsigned int i=0;i<N;i++) x[i] = 0.0;
  for(unsigned int i=0;i<N;i++) r[i] = b[i];

//  for( unsigned int i= 0; i< N; i++ )    // use defect correction
//    p[i] = r[i];

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
//    if( (k % 10) == 0 ) {
//      printf("    k,  |r|/|b| : %10d  %16.8E \n", k, cr );
  }
  return;
}
