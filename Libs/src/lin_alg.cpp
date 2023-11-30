#include "../include/lin_alg.hpp"

void linearBVP(double (*RHS)(const double& x), double y[], const double& xL, const double& xR, const double& yL, const double& yR, const int& nPoints) {
	y[0] = yL;
	y[nPoints - 1] = yR;
	const double dx = (xR - xL) / (nPoints - 1);

	double a[maxSize], b[maxSize], c[maxSize], r[maxSize];
	for (int i = 0; i < nPoints; i++) {
		a[i] =  1.0;
		b[i] = -2.0;
		c[i] =  1.0;
		r[i] = dx*dx * RHS(xL + i * dx);
	}
	r[1] -= y[0];
	r[nPoints - 2] -= y[nPoints - 1];
	// Required by tridiagonalSolver
	a[0] = a[1] = nan("");
	c[nPoints - 1] = c[nPoints - 2] = nan("");

	// Pass arrays + 1  (i.e.: arrays starting as index 1)
	//      nPoints - 2 (to stop at penultimate index)
	// Effectively, we're declaring the arrays 2 longer than we need them, then
	// scrapping the first and last element (the boundary conditions)
	tridiagonalSolver(a + 1, b + 1, c + 1, r + 1, y + 1, nPoints - 2);
}
