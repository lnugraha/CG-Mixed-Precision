make:
	g++ -std=c++17 -O3 -I./include ./src/Laplace.cpp ./src/MatrixOperations.cpp \
	./src/CGSolvers.cpp main.cpp -o cg_mixed_precision	
