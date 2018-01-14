double sdot(float* v0, float* v1, const int& N){
  float result = 0.0;
  for(int i=0;i<N;i++) result+=v0[i]*v1[i];
  return result;
}

void saxpyz(float* v0 , float* v1, const float& c , float* v2 , const int& N )
{
  for( unsigned int i= 0 ; i< N ; i++ )
    v0[i] = v1[i] + c * v2[i];

  return;
}

double ddot(double* v0, double* v1, const int& N){
  double result = 0.0;
  for(int i=0;i<N;i++) result+=v0[i]*v1[i];
  return result;
}

void daxpyz(double* v0 , double* v1, const double& c , double* v2 , const int& N )
{
  for( unsigned int i= 0 ; i< N ; i++ )
    v0[i] = v1[i] + c * v2[i];

  return;
}


