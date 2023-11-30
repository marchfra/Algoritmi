#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

void assignBoundaryConditions(double **M, double x[], double y[], const int& nx, const int& ny);

double boundaryCondition(const double& x, const double& y);

double SFunc(const double& x, const double& y);

int main() {
	const double xL = 0.0;
	const double xR = 1.0;
	const double yL = 0.0;
	const double yR = 1.0;
	const double tol = 1.0e-7;

	const int nx = 32, ny = 32;
	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	const double dx = (xR - xL) / nx;
	const double dy = (yR - yL) / ny;

	// Define matrices to store solution value at current and next iteration
	double **mOld;
	mOld = new double*[ny];
	mOld[0] = new double[nx * ny];
	for (int i = 1; i < ny; i++) {
		mOld[i] = mOld[i - 1] + nx;
	}
	double **mNew;
	mNew = new double*[ny];
	mNew[0] = new double[nx * ny];
	for (int i = 1; i < ny; i++) {
		mNew[i] = mNew[i - 1] + nx;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++)  x[i] = xL + i * dx;
	for (int j = 0; j < ny; j++)  y[j] = yL + j * dy;

	// Define source term and its value on the grid
	double **S;
	S = new double*[ny];
	S[0] = new double[nx * ny];
	for (int i = 1; i < ny; i++) {
		S[i] = S[i - 1] + nx;
	}
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			S[i][j] = SFunc(x[i], y[j]);
		}
	}

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++) {
		for (int j = jStart; j < jEnd; j++) {
			mOld[i][j] = 0.0;
		}
	}

	// Solve the equation
	double err = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition
		assignBoundaryConditions(mOld, x, y, nx, ny);

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				mNew[i][j] = 0.25 * (mOld[i + 1][j] + mOld[i - 1][j] + mOld[i][j + 1] + mOld[i][j - 1] - dx*dx * S[i][j]);
			}
		}

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				double d2Mdx2 = mOld[i + 1][j] - 2.0 * mOld[i][j] + mOld[i - 1][j];
				double d2Mdy2 = mOld[i][j + 1] - 2.0 * mOld[i][j] + mOld[i][j - 1];
				err += fabs(d2Mdx2 + d2Mdy2 - dx*dx * S[i][j]);
			}
		}

		// Assign new to old
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				mOld[i][j] = mNew[i][j];
			}
		}

		// Increment iteration counter
		numIter++;
	}

	cout << numIter << endl;

	// std::ofstream out;
	// out.open("../data/data.csv");
	// if (!out) exit(5);



	// out.close();

	delete[] mOld[0];
	delete[] mOld;
	delete[] mNew[0];
	delete[] mNew;

	return 0;
}

double boundaryCondition(const double& x, const double& y) {
	return exp(-M_PI * x) * sin(-M_PI * y) + 0.25 * SFunc(x, y) * (x*x + y*y);
}

double SFunc(const double& x, const double& y) {
	return 0.0;
}

void assignBoundaryConditions(double **M, double x[], double y[], const int& nx, const int& ny) {
	for (int i = 0,      j = 0;      j < ny; j++)  M[i][j] = boundaryCondition(x[i], x[j]);
	for (int i = nx - 1, j = 0;      j < ny; j++)  M[i][j] = boundaryCondition(x[i], x[j]);
	for (int i = 0,      j = 0;      i < nx; i++)  M[i][j] = boundaryCondition(x[i], x[j]);
	for (int i = 0,      j = ny - 1; i < nx; i++)  M[i][j] = boundaryCondition(x[i], x[j]);
}
