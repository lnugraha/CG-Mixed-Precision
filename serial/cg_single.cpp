/*
  Conjugate Gradient (CG) algorithm
*/

void  cg_single(float* x, float* b, float* p, int* site_ip, int* site_im, int* site_jp, int* site_jm, const int& N, const float& eps)
{

  float* r = new float[N];
  float* Ap = new float[N];        //  A.p
  float rr;                        //  (r,r)_k
  float tt;                        //  (r,r)_k+1
  float pAp;                       //  (p, A.p)
  float bb;                        //  (b, b)
  float alpha, beta;
  float cr;

// initialization
  for(unsigned int i=0;i<N;i++) x[i] = 0.0;  // Initial guess of x
  for(unsigned int i= 0;i<N;i++) r[i] = b[i]; // r = b - Ax_0 = b_0

  rr = sdot(r, r, N);	// Find rho
  bb = sdot(b, b, N);	// Find the |b|
  cr = sqrt(rr/bb);

  unsigned int k= 0;
  while (cr>eps){
    laplacian_times_vector_single(Ap, p, site_ip, site_im, site_jp, site_jm, N );      //   A.p
    pAp = sdot(p, Ap, N );                                               //   (p, Ap)
    alpha = (rr/pAp);
    saxpyz(x,  x,  alpha, p, N );         //  x = x + alpha * p
    saxpyz(r,  r, -alpha, Ap, N );        //  r = r - alpha * Ap
    tt = sdot(r, r, N);      //  tt = (r, r)
    beta = (tt/rr);
    saxpyz(p, r, beta, p, N);             //  p = r + beta * p
    k++;
    rr = tt;
    cr = sqrt(rr/bb);
  }
  return;
}
