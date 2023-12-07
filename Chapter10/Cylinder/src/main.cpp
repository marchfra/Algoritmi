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

void solveJacobi(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h);

void solveGaussSeidel(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h);

void solveSOR(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h, const double& omega);

int main() {
	const double xL = -1.0;
	const double xR =  1.0;
	const double yL = -1.0;
	const double yR =  1.0;

	const int nPoints = 128;
	const double tol = 1.0e-7;

	const double omega = 2.0 / (1.0 + M_PI / nPoints);

	const double h = (xR - xL) / (nPoints - 1);

	solveJacobi(nPoints, nPoints, xL, xR, yL, yR, tol, h);
	solveGaussSeidel(nPoints, nPoints, xL, xR, yL, yR, tol, h);
	solveSOR(nPoints, nPoints, xL, xR, yL, yR, tol, h, omega);

	return 0;
}

void solveJacobi(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h) {
	// ny = nRows
	// nx = nCols

	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrices to store solution value at current and next iteration and source matrix
	double **mOld;
	double **mNew;
	double **S;
	mOld = new double*[ny];
	mNew = new double*[ny];
	S = new double*[ny];
	mOld[0] = new double[ny * nx];
	mNew[0] = new double[ny * nx];
	S[0] = new double[ny * nx];
	for (int j = 1; j < ny; j++) {
		mOld[j] = mOld[j - 1] + nx;
		mNew[j] = mNew[j - 1] + nx;
		S[j] = S[j - 1] + nx;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++)  x[i] = xL + i * h;
	for (int j = 0; j < ny; j++)  y[j] = yL + j * h;

	// Assign source value on the grid
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
		// Assign boundary condition on mOld
		assignBoundaryConditions(mOld, x, y, nx, ny);

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				mNew[i][j] = 0.25 * (mOld[i + 1][j] + mOld[i - 1][j] + mOld[i][j + 1] + mOld[i][j - 1] - h*h * S[i][j]);
			}
		}

		// Assign boundary condition on mNew. Necessary for Laplacian error.
		assignBoundaryConditions(mNew, x, y, nx, ny);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// Laplacian error. Fine for Dirichlet conditions.
				double d2Mdx2 = mNew[i + 1][j] - 2.0 * mNew[i][j] + mNew[i - 1][j];
				double d2Mdy2 = mNew[i][j + 1] - 2.0 * mNew[i][j] + mNew[i][j - 1];
				err += fabs(d2Mdx2 + d2Mdy2 - h*h * S[i][j]);
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

	cout << "Number of iterations (Jacobi): " << numIter - 1 << endl;

	std::ofstream out;
	out.open("../data/data_jacobi.csv");
	if (!out) exit(5);

	out << "x,y,M,S" << endl;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			out << x[i] << "," << y[j] << "," << mOld[i][j] << "," << S[i][j] << endl;
		}
	}

	out.close();

	delete[] mOld[0];
	delete[] mOld;
	delete[] mNew[0];
	delete[] mNew;
	delete[] S[0];
	delete[] S;
}

void solveGaussSeidel(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h) {
	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrix to store solution values and source matrix
	double **M;
	double **S;
	M = new double*[ny];
	S = new double*[ny];
	M[0] = new double[ny * nx];
	S[0] = new double[ny * nx];
	for (int j = 1; j < ny; j++) {
		M[j] = M[j - 1] + nx;
		S[j] = S[j - 1] + nx;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++)  x[i] = xL + i * h;
	for (int j = 0; j < ny; j++)  y[j] = yL + j * h;

	// Assign source value on the grid
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			S[i][j] = SFunc(x[i], y[j]);
		}
	}

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++) {
		for (int j = jStart; j < jEnd; j++) {
			M[i][j] = 0.0;
		}
	}

	// Solve the equation
	double err = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition on M
		assignBoundaryConditions(M, x, y, nx, ny);

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				M[i][j] = 0.25 * (M[i + 1][j] + M[i - 1][j] + M[i][j + 1] + M[i][j - 1] - h*h * S[i][j]);
			}
		}

		// Assign boundary condition on M. Necessary for Laplacian error.
		assignBoundaryConditions(M, x, y, nx, ny);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// Laplacian error. Fine for Dirichlet conditions.
				double d2Mdx2 = M[i + 1][j] - 2.0 * M[i][j] + M[i - 1][j];
				double d2Mdy2 = M[i][j + 1] - 2.0 * M[i][j] + M[i][j - 1];
				err += fabs(d2Mdx2 + d2Mdy2 - h*h * S[i][j]);
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
			out << x[i] << "," << y[j] << "," << M[i][j] << "," << S[i][j] << endl;
		}
	}

	out.close();

	delete[] M[0];
	delete[] M;
	delete[] S[0];
	delete[] S;
}

void solveSOR(const int& nx, const int& ny, const double& xL, const double& xR, const double& yL, const double& yR, const double& tol, const double& h, const double& omega) {
	const int iStart = 1, iEnd = nx - 1;
	const int jStart = 1, jEnd = ny - 1;

	// Define matrix to store solution values and source matrix
	double **M;
	double **S;
	M = new double*[ny];
	S = new double*[ny];
	M[0] = new double[ny * nx];
	S[0] = new double[ny * nx];
	for (int j = 1; j < ny; j++) {
		M[j] = M[j - 1] + nx;
		S[j] = S[j - 1] + nx;
	}

	// Define grid points
	double x[nx], y[ny];
	for (int i = 0; i < nx; i++)  x[i] = xL + i * h;
	for (int j = 0; j < ny; j++)  y[j] = yL + j * h;

	// Assign source value on the grid
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			S[i][j] = SFunc(x[i], y[j]);
		}
	}

	// Initialize solution (initial guess)
	for (int i = iStart; i < iEnd; i++) {
		for (int j = jStart; j < jEnd; j++) {
			M[i][j] = 0.0;
		}
	}

	// Solve the equation
	double err = std::numeric_limits<double>::max();
	int numIter = 0;
	while (err > tol) {
		// Assign boundary condition on M
		assignBoundaryConditions(M, x, y, nx, ny);

		// Compute new value
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				M[i][j] = (1.0 - omega) * M[i][j] + 0.25 * omega * (M[i + 1][j] + M[i - 1][j] + M[i][j + 1] + M[i][j - 1] - h*h * S[i][j]);
			}
		}

		// Assign boundary condition on M. Necessary for Laplacian error.
		assignBoundaryConditions(M, x, y, nx, ny);

		// Compute error
		err = 0.0;
		for (int i = iStart; i < iEnd; i++) {
			for (int j = jStart; j < jEnd; j++) {
				// Laplacian error. Fine for Dirichlet conditions.
				double d2Mdx2 = M[i + 1][j] - 2.0 * M[i][j] + M[i - 1][j];
				double d2Mdy2 = M[i][j + 1] - 2.0 * M[i][j] + M[i][j - 1];
				err += fabs(d2Mdx2 + d2Mdy2 - h*h * S[i][j]);
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
			out << x[i] << "," << y[j] << "," << M[i][j] << "," << S[i][j] << endl;
		}
	}

	out.close();

	delete[] M[0];
	delete[] M;
	delete[] S[0];
	delete[] S;
}

double boundaryCondition(const double& x, const double& y) {
	const double a = 0.1;
	const double rho0 = 1.0;
	const double r = sqrt(x*x + y*y);
	if (r <= a) return -0.25 * rho0 * r*r;
	else		return -0.5  * rho0 * a*a * (log(r / a) + 0.5);
}

double SFunc(const double& x, const double& y) {
	const double a = 0.1;
	const double rho0 = 1.0;
	const double r = sqrt(x*x + y*y);
	if (r <= a) return -rho0;
	else		return 0.0;
}

void assignBoundaryConditions(double **M, double x[], double y[], const int& nx, const int& ny) {
	for (int i = 0,      j = 0;      j < ny; j++)  M[i][j] = boundaryCondition(x[i], y[j]);
	for (int i = nx - 1, j = 0;      j < ny; j++)  M[i][j] = boundaryCondition(x[i], y[j]);
	for (int i = 0,      j = 0;      i < nx; i++)  M[i][j] = boundaryCondition(x[i], y[j]);
	for (int i = 0,      j = ny - 1; i < nx; i++)  M[i][j] = boundaryCondition(x[i], y[j]);
}
