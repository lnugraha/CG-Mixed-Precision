/***************************
* Matrix vs Vector Toolbox *
* By: Leo Nugraha          *
* July 02, 2017	           *
***************************/
#ifndef _MVTOOLBOXH_
#define _MVTOOLBOXH_

// GENERAL MATRIX VECTOR MULTIPLICATION
void dgemv(double *sln, double *m, double *x, int size)
{
    int i,j; double sum;
    for(i=0;i<size;i++){
	sum = 0.0;
	for(j=0;j<size;j++) sum = sum + m[i*size+j]*x[j];
	sln[i] = sum;
    }
    return;
}

void sgemv(float *sln, float *m, float *x, int size)
{
   int i, j; float sum;
   for(i=0;i<size;i++){
	sum = 0.0;
	for(j=0;j<size;j++) sum = sum + m[i*size+j]*x[j];
	sln[i] = sum;
   }
   return;
}

// VECTOR CONVERSION
void svec_add_dvec(float *in, double *out, int size)
{
    for(int i=0;i<size;i++) out[i] += (double)in[i];
    return;
}


// VECTOR SCALING OPERATIONS
void vec_add_scaled_vec(double *sln, double *a, double c, double *b, int size)
{
    int i;
    for(i=0; i<size; i++) sln[i] = a[i]+c*b[i];
    return;
}

void vec_add_scaled_vec_sgl(float *sln, float *a, float c, float *b, int size)
{
    int i;
    for(i=0;i<size;i++) sln[i] = a[i]+c*b[i];
    return;
}

// Handle dot products between 2 vectors
// out_vector = dot_product(vector_01, vector_02, length)
double ddot(double *a, double *b, int size)
{
    double sum = 0.0;
    for(int i=0; i<size; i++) sum = sum+a[i]*b[i];
    return sum;
}

float sdot(float *a, float *b, int size)
{
    float sum = 0.0;
    for(int i=0;i<size;i++) sum=sum+a[i]*b[i];
    return sum;
}

// VECTOR CLONE OPERATIONS
void vector_dcopy(double *copy, double *orig, int size)
{
    for(int i=0;i<size;i++)
	copy[i] = orig[i];
}

void vector_scopy(float *copy, float *orig, int size)
{
    for(int i=0;i<size;i++) copy[i] = orig[i];
}


#endif
