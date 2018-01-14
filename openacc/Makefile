main : *.cpp
	pgc++ -Minfo=accel -ta=multicore -acc -fast -O3 -o dbl_gpu gpu_solver.cpp -Msafeptr
clean :
	rm dbl_omp_jan12
