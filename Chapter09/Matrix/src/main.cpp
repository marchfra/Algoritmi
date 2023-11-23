#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#define NROW 4
#define NCOL 4

void matrixVectorMult(double **M, double v[], const int& nRowsM, const int& nColsM, const int& nRowsV, double result[]);

int main() {
	// std::ofstream out;
	// out.open("data/data.csv");
	// if (!out) exit(5);
	
	double **M;
	M = new double*[NROW];
	M[0] = new double[NROW*NCOL];
	for (int i = 1; i < NROW; i++) {
		M[i] = M[i - 1] + NCOL;
	}

	M[0][0] = 1;	M[0][1] = 3;	M[0][2] = 2;	M[0][3] = -4;
	M[1][0] = 7;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 1;
	M[2][0] = 0;	M[2][1] = -1;	M[2][2] = 2;	M[2][3] = 2;
	M[3][0] = 6;	M[3][1] = 3;	M[3][2] = 0;	M[3][3] = 1;

	cout.setf(std::ios::fixed);
	cout.precision(4);

	for (int i = 0; i < NROW; i++) {
		for (int j = 0; j < NCOL; j++) {
			cout << std::setw(8) << M[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	double b[NCOL];

	b[0] = 1;
	b[1] = 0;
	b[2] = 3;
	b[3] = 2;

	for (int i = 0; i < NCOL; i++) {
		cout << std::setw(8) << b[i] << endl;
	}

	// cout << endl << "Multiplying M * b" << endl;
	// for (int i = 0; i < NROW; i++) {
	// 	double s = 0.0;
	// 	for (int j = 0; j < NCOL; j++) {
	// 		s += M[i][j] * b[j];
	// 		// cout << M[i][j] << " * " << b[j] << " = ";
	// 		// cout << s << "; ";
	// 	}
	// 	cout << std::setw(8) << s << endl;
	// }

	cout << endl;
	double result[NCOL];

	matrixVectorMult(M, b, NROW, NCOL, NCOL, result);

	for (int i = 0; i < NCOL; i++) {
		cout << std::setw(8) << result[i] << endl;
	}

	// cout << "Hello World!\n";

	delete[] M[0];
	delete[] M;

	// out.close();
	return 0;
}

void matrixVectorMult(double **M, double *v, const int& nRowsM, const int& nColsM, const int& nRowsV, double result[]) {
	if (nRowsM != nRowsV) {
		throw std::invalid_argument("matrixVectorMult(): nRowsMatrix must be nRowsVector");
	}

	for (int i = 0; i < NROW; i++) {
		double s = 0.0;
		for (int j = 0; j < NCOL; j++) {
			s += M[i][j] * v[j];
		}
		result[i] = s;
		// cout << std::setw(8) << s << endl;
	}
}


double f(const double& x);

double f(double x)
