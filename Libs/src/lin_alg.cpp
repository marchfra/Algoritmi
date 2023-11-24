#include <iostream>
#include <iomanip>
#include <string>

#include "../include/lin_alg.hpp"
#include "../include/swap.hpp"

void printMatrix(const double **M, const int& n, const int width) {
	for (int i = 0; i < n; i++) {
		std::cout << std::string(width + 1, '-');
	}
	std::cout << std::endl;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << std::setw(width) << M[i][j] << " ";
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < n; i++) {
		std::cout << std::string(width + 1, '-');
	}
	std::cout << std::endl;
}

void printVector(const double v[], const int& nRows) {
	std::cout << "{";
	for (int i = 0; i < nRows; i++) {
		std::cout << v[i];
		if (i != nRows - 1) std::cout << ", ";
	}
	std::cout << "}" << std::endl;
}

void swapRowsLinSystem(double **M, double v[], const int& nRows, const int& i, const int& j) {
	for (int k = 0; k < nRows; k++) {
		swap(M[i][k], M[j][k]);
	}

	swap(v[i], v[j]);
}

void gaussElim(double **M, double v[], const int& nRows) {
	// Reduce M to upper triangular form
	for (int k = 0; k < nRows - 1; k++) {
		partialPivoting(M, v, nRows, k);
		for (int i = k + 1; i < nRows; i++) {
			double g = M[i][k] / M[k][k];
			for (int j = k + 1; j < nRows; j++) M[i][j] -= g * M[k][j];
			M[i][k] = 0.0;
			v[i] -= g * v[k];
		}
	}
}

void partialPivoting(double **M, double v[], const int& nRows, const int& current_row) {
	double max = fabs(M[current_row][current_row]), maxI = current_row;
	for (int i = current_row; i < nRows; i++) {	// Loop over rows under current row
		double g = fabs(M[i][current_row]);
		if (g > max) {
			max = g;
			maxI = i;
		}
	}
	if (maxI != current_row) {
		swapRowsLinSystem(M, v, nRows, maxI, current_row);
	}
}

void solveLinSystem(double **M, double v[], double x[], const int& nEqs) {
	// Reduce the system to upper triangular form
	gaussElim(M, v, nEqs);

	// Solve the system
	for (int i = nEqs - 1; i >= 0; i--) {
		double temp = v[i];
		for (int j = nEqs - 1; j > i; j--) temp -= x[j] * M[i][j];
		x[i] = temp / M[i][i];
	}
}
