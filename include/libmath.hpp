#ifndef LIBMATH_HPP_
#define LIBMATH_HPP_

class single{
	public:
		float sdot(float *v0, float *v1, const int &N);
		void saxpyz(float *v0, float *v1, float *v2, const float &c, const int &N);
};

class dual{
	public:
		double ddot(double *v0, double *v1, const int &N);
		void daxpyz(double *v0, double *v1, double *v2, const double &c, const int &N);
};

class mutual{
	public:
		void vec_add(double *vec_dbl, const float *vec_sgl, const int &N);
};

#endif
