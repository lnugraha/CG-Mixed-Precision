void vector_addition_dbl_sgl( double* x_dbl , const float* x_sgl , const int& N )
{
  for( unsigned int i= 0; i< N; i++)
    x_dbl[i] += (double) x_sgl[i];

  return;
}
