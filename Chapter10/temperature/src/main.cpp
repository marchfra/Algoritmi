#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
// #include
// "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

void assignBoundaryConditions(double** M, double x[], double y[], const int& nx,
                              const int& ny, const double& h);

double SFunc(const double& x, const double& y);

void solveJacobi(const int& nx, const int& ny, const double& xL,
                 const double& xR, const double& yL, const double& yR,
                 const double& tol, const double& h);

void solveGaussSeidel(const int& nx, const int& ny, const double& xL,
                      const double& xR, const double& yL, const double& yR,
                      const double& tol, const double& h);

void solveSOR(const int& nx, const int& ny, const double& xL, const double& xR,
              const double& yL, const double& yR, const double& tol,
              const double& h, const double& omega);

int main() {
	const double xL = 0.0;
	const double xR = 2.0;
	const double yL = 0.0;
	const double yR = 1.0;

	const int nx     = 129;
	const int ny     = 65;
	const double tol = 1.0e-7;

	const double omega = 1.0;  // 2.0 / (1.0 + M_PI / nPoints);

	const double hx = (xR - xL) / (nx - 1);
	const double hy = (yR - yL) / (ny - 1);
	const double h =
		(hx == hy ? hx
	              : throw std::invalid_argument("hx must be equal to hy."));

	solveJacobi(nx, ny, xL, xR, yL, yR, tol, h);
	solveGaussSeidel(nx, ny, xL, xR, yL, yR, tol, h);
	solveSOR(nx, ny, xL, xR, yL, yR, tol, h, omega);

	return 0;
}

void solveJacobi(const int& nx, const int& ny, const double& xL,
                 const double& xR, const double& yL, const double& yR,
                 const double& tol, const double& h) {
	// nx = nRows
	// ny = nCols

	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrices to store solution value at current and next iteration and
	// source matrix
	double** mOld;
	double** M;
	double** S;
	mOld    = new double*[nx];
	M       = new double*[nx];
	S       = new double*[nx];
	mOld[0] = new double[nx * ny];
	M[0]    = new double[nx * ny];
	S[0]    = new double[nx * ny];
	for (int j = 1; j < nx; j++) {
		mOld[j] = mOld[j - 1] + ny;
		M[j]    = M[j - 1] + ny;
		S[j]    = S[j - 1] + ny;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++) x[i] = xL + i * h;
	for (int j = 0; j < ny; j++) y[j] = yL + j * h;

	// Assign source value on the grid
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) S[i][j] = SFunc(x[i], y[j]);

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++)
		for (int j = jStart; j < jEnd; j++) mOld[i][j] = 0.0;

	// Solve the equation
	double err  = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition on mOld
		assignBoundaryConditions(mOld, x, y, nx, ny, h);

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				M[i][j] =
					0.25 * (mOld[i + 1][j] + mOld[i - 1][j] + mOld[i][j + 1] +
				            mOld[i][j - 1] - h * h * S[i][j]);
			}
		}

		// Assign boundary condition on M. Necessary for Laplacian error.
		assignBoundaryConditions(M, x, y, nx, ny, h);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// Laplacian error. Fine for Dirichlet conditions.
				// double d2Mdx2 =
				// 	M[i + 1][j] - 2.0 * M[i][j] + M[i - 1][j];
				// double d2Mdy2 =
				// 	M[i][j + 1] - 2.0 * M[i][j] + M[i][j - 1];
				// err += fabs(d2Mdx2 + d2Mdy2 - h * h * S[i][j]);

				// "Convergence error". Better for Neumann conditions.
				err += fabs(M[i][j] - mOld[i][j]) * h * h;
			}
		}

		// Assign new to old
		for (int i = iStart; i < iEnd; i++)
			for (int j = jStart; j < jEnd; j++) mOld[i][j] = M[i][j];

		// Increment iteration counter
		numIter++;
	}

	cout << "Number of iterations (Jacobi): " << numIter - 1 << endl;

	std::ofstream out;
	out.open("../data/data_jacobi.csv");
	if (!out) exit(5);

	out << "x,y,M,S" << endl;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << x[i] << "," << y[j] << "," << mOld[i][j] << "," << S[i][j]
				<< endl;
		}
	}

	out.close();

	delete[] mOld[0];
	delete[] mOld;
	delete[] M[0];
	delete[] M;
	delete[] S[0];
	delete[] S;
}

void solveGaussSeidel(const int& nx, const int& ny, const double& xL,
                      const double& xR, const double& yL, const double& yR,
                      const double& tol, const double& h) {
	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrix to store solution values and source matrix
	double** mOld;
	double** M;
	double** S;
	mOld    = new double*[nx];
	M       = new double*[nx];
	S       = new double*[nx];
	mOld[0] = new double[nx * ny];
	M[0]    = new double[nx * ny];
	S[0]    = new double[nx * ny];
	for (int j = 1; j < nx; j++) {
		mOld[j] = mOld[j - 1] + ny;
		M[j]    = M[j - 1] + ny;
		S[j]    = S[j - 1] + ny;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++) x[i] = xL + i * h;
	for (int j = 0; j < ny; j++) y[j] = yL + j * h;

	// Assign source value on the grid
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) S[i][j] = SFunc(x[i], y[j]);

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++)
		for (int j = jStart; j < jEnd; j++) M[i][j] = 0.0;

	// Solve the equation
	double err  = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition on M
		assignBoundaryConditions(M, x, y, nx, ny, h);

		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++) mOld[i][j] = M[i][j];

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				M[i][j] = 0.25 * (M[i + 1][j] + M[i - 1][j] + M[i][j + 1] +
				                  M[i][j - 1] - h * h * S[i][j]);
			}
		}

		// Assign boundary condition on M. Necessary for Laplacian error.
		assignBoundaryConditions(M, x, y, nx, ny, h);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// "Convergence error". Better for Neumann conditions.
				err += fabs(M[i][j] - mOld[i][j]) * h * h;
			}
		}

		// Increment iteration counter
		numIter++;
	}

	cout << "Number of iterations (Gauss - Seidel): " << numIter - 1 << endl;

	std::ofstream out;
	out.open("../data/data_gauss.csv");
	if (!out) exit(5);

	out << "x,y,M,S" << endl;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << x[i] << "," << y[j] << "," << M[i][j] << "," << S[i][j]
				<< endl;
		}
	}

	out.close();

	delete[] mOld[0];
	delete[] mOld;
	delete[] M[0];
	delete[] M;
	delete[] S[0];
	delete[] S;
}

void solveSOR(const int& nx, const int& ny, const double& xL, const double& xR,
              const double& yL, const double& yR, const double& tol,
              const double& h, const double& omega) {
	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrix to store solution values and source matrix
	double** mOld;
	double** M;
	double** S;
	mOld    = new double*[nx];
	M       = new double*[nx];
	S       = new double*[nx];
	mOld[0] = new double[nx * ny];
	M[0]    = new double[nx * ny];
	S[0]    = new double[nx * ny];
	for (int j = 1; j < nx; j++) {
		mOld[j] = mOld[j - 1] + ny;
		M[j]    = M[j - 1] + ny;
		S[j]    = S[j - 1] + ny;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++) x[i] = xL + i * h;
	for (int j = 0; j < ny; j++) y[j] = yL + j * h;

	// Assign source value on the grid
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) S[i][j] = SFunc(x[i], y[j]);

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++)
		for (int j = jStart; j < jEnd; j++) M[i][j] = 0.0;

	// Solve the equation
	double err  = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition on M
		// assignBoundaryConditions(M, x, y, nx, ny, h);

		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++) mOld[i][j] = M[i][j];

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				M[i][j] = (1.0 - omega) * M[i][j] +
				          0.25 * omega *
				              (M[i + 1][j] + M[i - 1][j] + M[i][j + 1] +
				               M[i][j - 1] - h * h * S[i][j]);
			}
		}

		// Assign boundary condition on M. Necessary for Laplacian error.
		assignBoundaryConditions(M, x, y, nx, ny, h);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// "Convergence error". Better for Neumann conditions.
				err += fabs(M[i][j] - mOld[i][j]) * h * h;
			}
		}

		// Increment iteration counter
		numIter++;
	}

	cout << "Number of iterations (SOR): " << numIter - 1 << endl;

	std::ofstream out;
	out.open("../data/data_sor.csv");
	if (!out) exit(5);

	out << "x,y,M,S" << endl;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << x[i] << "," << y[j] << "," << M[i][j] << "," << S[i][j]
				<< endl;
		}
	}

	out.close();

	delete[] mOld[0];
	delete[] mOld;
	delete[] M[0];
	delete[] M;
	delete[] S[0];
	delete[] S;
}

double SFunc(const double& x, const double& y) {
	return 0.0;
}

void assignBoundaryConditions(double** M, double x[], double y[], const int& nx,
                              const int& ny, const double& h) {
	for (int i = 0, j = 0; j < ny; j++)  // b.c. at x = xL
		M[i][j] = M[i + 1][j] - 0.0 * h;
	for (int i = nx - 1, j = 0; j < ny; j++)  // b.c. at x = xR
		M[i][j] = M[i - 1][j] + 3.0 * h;
	for (int i = 0, j = 0; i < nx; i++)  // b.c. at y = yL
		M[i][j] = 0.0;
	for (int i = 0, j = ny - 1; i < nx; i++)  // b.c. at y = yR
		M[i][j] = 2.0 - x[i];
}
