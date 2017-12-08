/******************************************************************
* CG Iteration Method - Basic Mode				  *
* Input: A Symmetric, Positive, Definite SQUARE matrix		  *
* Outputs: Number of Iterations to achieve convergence & Errors   *
* (Note) You can also print the solution			  *
* By: Leo Nugraha						  *
* Date: Nov 08, 2017						  *
******************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mvtoolbox.h"
#include "cg_double_solve.cpp"

#include "mmio.h"
#include "mmio.c"

clock_t start,end;	// CPU Timing - from time.h
double cpu_time_used;	// Compute the time difference

int main(int argc, char *argv[])
{
    int ret_code, size, n, M, N, nz, i, *I, *J;
    MM_typecode matcode;
    FILE *f;
    double *value;
    FILE *inputfile;

    if (argc < 2) {
	fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
	exit(1);
    }
    else {
        if ((f = fopen(argv[1], "r")) == NULL)
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)	exit(1);

    /* reseve memory for matrices */
    I = (int*)malloc(nz*sizeof(int));
    J = (int*)malloc(nz*sizeof(int));
    value = (double*)malloc(nz*sizeof(double));

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &value[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);
	size = M;
    printf("Matrix Size: %d \n", size);
    printf("Non-Zero Elements: %d \n", nz);

    // Define the solution b = [1 1 1 ... 1]'
    double *b;
    b = (double*)malloc(size*sizeof(double));
    for (i = 0; i < size; i++)	b[i] = 1.0;

    // Define the input A = (size X size), where A is an SPD matrix
    // Note: To ensure the convergency, A must be a diagonal dominant matrix (the matrix diagonal has the highest value)
    double *A;
    A = (double*)malloc(size*size*sizeof(double*));
    for(i = 0; i < size*size; i++)	A[i] = 0.0;

    for(int a=0;a<nz;a++){
    	int n = I[a]*size+J[a];
    	A[n] = value[a];
    }
    for(int a=0;a<nz;a++){
    	int m = J[a]*size+I[a];
    	A[m] = value[a];
    }

    // Define the output x = [... ... ... ...]', which will be determined through iterations
    // Our initial guess is all zeroes, and will be updated for each iteration
    double *x;
    x = (double*)malloc(size*sizeof(double));
    for(i = 0; i < size; i++)	x[i] = 0.0;

    double eps   = 1.0e-6;	// Error Tolerance
    int maxiter  = size*size;	// Maximum Number of Iteration
    int cnt      = 0;		// Number of iteration used

    start = clock();
    cg_double(size, A, b, eps, maxiter, &cnt, x);
    end = clock();

    printf("Computed %d Iterations\n",cnt);
    /* Compute the Error */
    double sum = 0.0;
    for(i=0; i<size; i++){
	double d = x[i] - 1.0;
	sum += (d >= 0.0) ? d : -d;
    }
    printf("Error : %.3e\n",sum);

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf( "CPU time used %16.8E seconds \n",cpu_time_used);

    // See your solution
//     for (i = 0; i < size; i++)
//    	printf("%.6lf \n", x[i]);
   free(I); free(J); free(value);
   free(b); free(x); free(A);
   return 0;
}
