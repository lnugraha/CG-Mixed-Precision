#ifndef LAPLACIAN_HPP_
#define LAPLACIAN_HPP_


class laplace{
	public:
	void siteindex(int *site_ip, int *site_im, int *site_jp, int *site_jm, const int &N);

	void laplacian_times_vector( double* Ax, const double* x, int* site_ip, int* site_im, 
                             int* site_jp, int* site_jm, const int& N );

	void laplacian_times_vector_single( float* Ax, const float* x, int* site_ip, int* site_im, 
                                    int* site_jp, int* site_jm, const int& N );
};

#endif
