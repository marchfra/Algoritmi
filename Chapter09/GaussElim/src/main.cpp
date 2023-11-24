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

// void gaussElim(double **M, double v[], double x[], const int& n);

// void printM(double **M, const int& n, const int width=8);

// void printV(double v[], const int& n);

// void swapRows(double **M, double v[], const int& n, const int& i, const int& j);

int main() {
	// static const int N = 4;
	// double **M;
	// M = new double*[N];
	// M[0] = new double[N*N];
	// for (int i = 1; i < N; i++) {
	// 	M[i] = M[i - 1] + N;
	// }
	// M[0][0] = 1;	M[0][1] = 2;	M[0][2] = 1;	M[0][3] = -1;
	// M[1][0] = 3;	M[1][1] = 2;	M[1][2] = 4;	M[1][3] = 4;
	// M[2][0] = 4;	M[2][1] = 4;	M[2][2] = 3;	M[2][3] = 4;
	// M[3][0] = 2;	M[3][1] = 0;	M[3][2] = 1;	M[3][3] = 5;

	// double v[N] = {5, 16, 22, 15};

	// printM(M, N, 4);
	// printV(v, N);

	// swapRows(M, v, N, 0, 2);

	// printM(M, N, 4);
	// printV(v, N);

	// delete[] M[0];
	// delete[] M;


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
	cout << "Solution of the first system:  x = ";
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
	return 0;
}

// void gaussElim(double **M, double v[], double x[], const int& n) {
// 	// Reduce M to upper triangular form
// 	for (int k = 0; k < n - 1; k++) {
// 		// Partial pivoting
// 		double max = fabs(M[k][k]), maxI = k;
// 		for (int i = k; i < n; i++) {
// 			if (fabs(M[i][k]) > max) {
// 				max = fabs(M[i][k]);
// 				maxI = i;
// 			}
// 		}
// 		if (maxI != k) {
// 			swapRows(M, v, n, maxI, k);
// 		}
// 		for (int i = k + 1; i < n; i++) {
// 			double g = M[i][k] / M[k][k];
// 			for (int j = k + 1; j < n; j++) M[i][j] -= g * M[k][j];
// 			M[i][k] = 0.0;
// 			v[i] -= g * v[k];
// 		}
// 	}

// 	// Solve the system
// 	for (int i = n - 1; i >= 0; i--) {
// 		double temp = v[i];
// 		for (int j = n - 1; j > i; j--) temp -= x[j] * M[i][j];
// 		x[i] = temp / M[i][i];
// 	}
// }


// void printM(double **M, const int& n, const int width) {
// 	for (int i = 0; i < n; i++) {
// 		cout << std::string(width + 1, '-');
// 	}
// 	cout << endl;

// 	for (int i = 0; i < n; i++) {
// 		for (int j = 0; j < n; j++) {
// 			cout << std::setw(width) << M[i][j] << " ";
// 		}
// 		cout << endl;
// 	}

// 	for (int i = 0; i < n; i++) {
// 		cout << std::string(width + 1, '-');
// 	}
// 	cout << endl;
// }

// void printV(double v[], const int& n) {
// 	cout << "{";
// 	for (int i = 0; i < n; i++) {
// 		cout << v[i];
// 		if (i != n - 1) cout << ", ";
// 	}
// 	cout << "}" << endl;
// }

// void swapRows(double **M, double v[], const int& n, const int& i, const int& j) {
// 	for (int k = 0; k < n; k++) {
// 		swap(M[i][k], M[j][k]);
// 	}

// 	swap(v[i], v[j]);
// }
