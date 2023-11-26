#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

int main() {
	// First system
	static const int N1 = 3;

	double **M1;
	M1 = new double*[N1];
	M1[0] = new double[N1*N1];
	for (int i = 1; i < N1; i++) {
		M1[i] = M1[i - 1] + N1;
	}

	M1[0][0] = 2;	M1[0][1] = -1;	M1[0][2] = 0;
	M1[1][0] = -1;	M1[1][1] = 2;	M1[1][2] = -1;
	M1[2][0] = 0;	M1[2][1] = -1;	M1[2][2] = 2;

	double v1[N1] = {1, 2, 1};
	double x1[N1];

	solveLinSystem(M1, v1, x1, N1);
	// printMatrix(M1, N1);
	cout << "Solution of the first system: x = ";
	printVector(x1, N1);

	delete[] M1[0];
	delete[] M1;


	// Second system
	static const int N2 = 4;

	double **M2;
	M2 = new double*[N2];
	M2[0] = new double[N2*N2];
	for (int i = 1; i < N2; i++) {
		M2[i] = M2[i - 1] + N2;
	}

	M2[0][0] = 1;	M2[0][1] = 2;	M2[0][2] = 1;	M2[0][3] = -1;
	M2[1][0] = 3;	M2[1][1] = 2;	M2[1][2] = 4;	M2[1][3] = 4;
	M2[2][0] = 4;	M2[2][1] = 4;	M2[2][2] = 3;	M2[2][3] = 4;
	M2[3][0] = 2;	M2[3][1] = 0;	M2[3][2] = 1;	M2[3][3] = 5;

	double v2[N2] = {5, 16, 22, 15};
	double x2[N2];

	solveLinSystem(M2, v2, x2, N2);
	// printMatrix(M2, N2);
	cout << "Solution of the second system: x = ";
	printVector(x2, N2);

	delete[] M2[0];
	delete[] M2;


	// Third system
	static const int N3 = 4;

	double **M3;
	M3 = new double*[N3];
	M3[0] = new double[N3*N3];
	for (int i = 1; i < N3; i++) {
		M3[i] = M3[i - 1] + N3;
	}

	M3[0][0] = 1;	M3[0][1] = 2;	M3[0][2] = 1;	M3[0][3] = -1;
	M3[1][0] = 3;	M3[1][1] = 6;	M3[1][2] = 4;	M3[1][3] = 4;
	M3[2][0] = 4;	M3[2][1] = 4;	M3[2][2] = 3;	M3[2][3] = 4;
	M3[3][0] = 2;	M3[3][1] = 0;	M3[3][2] = 1;	M3[3][3] = 5;

	double v3[N3] = {5, 16, 22, 15};
	double x3[N3];

	solveLinSystem(M3, v3, x3, N3);
	cout << "Solution of the third system:  x = ";
	printVector(x3, N3);

	delete[] M3[0];
	delete[] M3;


	// Tridiagonal system
	const int nEq = 5;
	double a[nEq] = {nan(""), 1, 1, 1, 1};
	double b[nEq] = {2, 2, 2, 2, 2};
	double c[nEq] = {1, 1, 1, 1, nan("")};
	double v[nEq] = {1, 0, 3, 1, 0};
	double x[nEq];

	tridiagonalSolver(a, b, c, v, x, nEq);
	cout << "Solution of the tridiagonal system: x = ";
	printVector(x, nEq);

	return 0;
}
