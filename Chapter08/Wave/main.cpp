#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

void RHS(double, double[], double[]);

double Residual(double);

int main() {
	// Shooting
	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);

	out << "x,phi,k" << endl;

	const double xa = 0.0, xb = 1.0;
	const int nStep = 100;
	const double dx = (xb - xa) / nStep;
	double x = xa;

	const double phi0 = 0.0;
	const double dphi0 = 1.0;
	const double kmin = 0.0, kmax = 5.0;
	
	for (double k = kmin; k <= kmax; k++) {
		double y[] = {phi0, dphi0, k};
		const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

		x = xa;
		out << x << "," << y[0] << "," << y[2] << endl;
		for (int i = 0; i < nStep; i++) {
			RK4Step(x, y, RHS, dx, nEq);
			x += dx;

			out << x << "," << y[0] << "," << y[2] << endl;
		}
	}

	out.close();

	// Find the roots of Residual(k) for k in [1, 20]
	out.open("final.csv");
	if (!out) exit(1);

	const double ka =  1.0;
	const double kb = 20.0;
	const double tol = 1.0e-8;
	double k_roots[64];
	int nRoots;
	find_roots(Residual, ka, kb, tol, k_roots, nRoots);

	out << "x,phi,k" << endl;
	for (int i = 0; i < nRoots; i++) {
		cout << "root #" << i << ": " << k_roots[i] / M_PI << "π" << endl;

		double y_final[] = {phi0, dphi0, k_roots[i]};
		const int nEq = static_cast<int>(sizeof(y_final)) / static_cast<int>(sizeof(y_final[0]));

		x = xa;
		out << x << "," << y_final[0] << "," << y_final[2] << endl;
		for (int i = 0; i < nStep; i++) {
			RK4Step(x, y_final, RHS, dx, nEq);
			x += dx;

			out << x << "," << y_final[0] << "," << y_final[2] << endl;
		}
	}

	out.close();

	return 0;
}

void RHS(double x, double Y[], double R[]) {
	double phi = Y[0];
	double dphi = Y[1];
	double k = Y[2];

	R[0] = dphi;
	R[1] = -k*k * phi;
	R[2] = 0;
}

double Residual(double k) {
	const double xa = 0.0, xb = 1.0;
	const int nStep = 100;
	const double dx = (xb - xa) / nStep;
	double x = xa;

	const double phi0 = 0.0;
	const double dphi0 = 1.0;

	double y[] = {phi0, dphi0, k};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx, nEq);
		x += dx;
	}

	return y[0];
}
