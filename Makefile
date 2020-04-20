mixed_precision:
	g++ -std=c++17 -O3 -I./include ./src/Laplace.cpp ./src/MatrixOperations.cpp \
	./src/CGSolvers.cpp cg_mixed.cpp -o mixed_prec	

double_precision:
	g++ -std=c++17 -O3 -I./include ./src/Laplace.cpp ./src/MatrixOperations.cpp \
	./src/CGSolvers.cpp cg_double.cpp -o double_prec
