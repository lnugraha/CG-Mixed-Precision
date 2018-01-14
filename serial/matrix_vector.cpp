void laplacian_times_vector( double* Ax, const double* x, int* site_ip, int* site_im, 
                             int* site_jp, int* site_jm, const int& N )
{
  double mass = 0.1;

  for( unsigned int i= 0 ; i< N ; i++ ) 
    Ax[i] = (mass*mass-4.0)*x[i] + x[site_ip[i]] + x[site_im[i]] + x[site_jp[i]] + x[site_jm[i]];
  
  return;
}

void laplacian_times_vector_single( float* Ax, const float* x, int* site_ip, int* site_im, 
                                    int* site_jp, int* site_jm, const int& N )
{
  float mass = 0.1;

  for( unsigned int i= 0 ; i< N ; i++ ) 
    Ax[i] = (mass*mass-4.0)*x[i] + x[site_ip[i]] + x[site_im[i]] + x[site_jp[i]] + x[site_jm[i]];
  
  return;
}
