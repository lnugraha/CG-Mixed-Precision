#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mvtoolbox.h"

// void vector_dcopy(double *, double *, int);
void cg_double(int size, double *A, double *b, double epsilon, int maxiter, int *numit, double *x){
    double alpha, beta, rdr, new_rdr, dMd, bdb;
    double *r_ptr = (double*)malloc(size*sizeof(double));
    double *d_ptr = (double*)malloc(size*sizeof(double));
    double *e_ptr = (double*)malloc(size*sizeof(double));

    dgemv(e_ptr, A, x, size);
    vec_add_scaled_vec(d_ptr, b, -1.0, e_ptr, size);
    vector_dcopy(r_ptr, d_ptr, size);
    rdr = ddot(r_ptr, r_ptr, size);
    bdb = ddot(b, b, size);//NEW
    double cr = sqrt(rdr/bdb);//NEW

    while(cr > epsilon){
	dgemv(e_ptr, A, d_ptr, size);
	dMd = ddot(d_ptr, e_ptr, size);
	alpha = rdr/dMd;

	vec_add_scaled_vec(x,x,alpha,d_ptr,size);
	vec_add_scaled_vec(r_ptr,r_ptr,-alpha,e_ptr,size);

	new_rdr = ddot(r_ptr,r_ptr,size);
	beta = new_rdr/rdr;

	vec_add_scaled_vec(d_ptr,r_ptr,beta,d_ptr,size);
	rdr = new_rdr;

	bdb = ddot(b,b,size);
	cr = sqrt(rdr/bdb);//NEW

	*numit = *numit + 1;
	int k = *numit;
	if (k>maxiter){
	   printf("FAIL TO CONVERGE! \n");
	   exit(1);
	}
    }
free(r_ptr); free(d_ptr); free(e_ptr);
}

/*
void vector_dcopy(double *copy, double *orig, int size){
	for(int i=0;i<size;i++) copy[i]=orig[i];

} */
