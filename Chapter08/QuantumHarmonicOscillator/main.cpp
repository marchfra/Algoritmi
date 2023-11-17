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

double V(double);

void RHS(double, double[], double[]);

void ForwardInt();

void BackwardInt();

void MatchingPointShooting();

double Res(double);

void ResPlot();

double ResZero();

void Final();

void reset_y(double[], const int&, const double&, const double&, const double&);

int main() {
	ForwardInt();
	BackwardInt();
	MatchingPointShooting();
	ResPlot();
	Final();

	return 0;
}

double V(double x) {
	return 0.5 * x*x;
}

void RHS(double x, double Y[], double R[]) {
	double psi = Y[0];
	double dpsi = Y[1];
	double E = Y[2];

	R[0] = dpsi;					// dpsi/dx = s
	R[1] = -2.0 * (E - V(x)) * psi;	// d2psi/dx2 = -2(V - E) psi
	R[2] = 0.0;						// dE/dx = 0;
}

void ForwardInt() {
	std::ofstream out;
	out.open("forward.csv");
	if (!out) exit(1);

	out << "x,psi" << endl;

	const double xa = -10.0, xb = 10.0;
	const int nStep = 800;
	const double dx = (xb - xa) / nStep;
	double x = xa;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = -xa * psi0;
	const double E0 = 0.5;
	double y[] = {psi0, dpsi0, E0};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx, nEq);
		x += dx;

		out << x << "," << y[0] << endl;
	}

	out.close();
}

void BackwardInt() {
	std::ofstream out;
	out.open("backward.csv");
	if (!out) exit(1);

	out << "x,psi" << endl;

	const double xa = 10.0, xb = -10.0;
	const int nStep = 800;
	const double dx = (xb - xa) / nStep;
	double x = xa;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = -xa * psi0;
	const double E0 = 0.5;
	double y[] = {psi0, dpsi0, E0};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx, nEq);
		x += dx;

		out << x << "," << y[0] << endl;
	}

	out.close();
}

void MatchingPointShooting() {
	std::ofstream out;
	out.open("matching_shoot.csv");
	if (!out) exit(1);

	out << "x,psi,E" << endl;

	const double xa = -10.0, xb = 10.0;
	const double xm = 0.5 * (xa + xb) + 1e-1;
	const int nStep = 800;
	const double dx1 = (xm - xa) / nStep;
	const double dx2 = (xm - xb) / nStep;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = -xa * psi0;
	const double psi1 = exp(-0.5 * xb*xb);
	const double dpsi1 = -xb * psi1;
	double E = 0.0;
	double y[] = {psi0, dpsi0, E};
	const int nEq = static_cast<int>(sizeof(y) / sizeof(y[0]));

	while (E < 5.0) {
		// Forward integration
		reset_y(y, nEq, psi0, dpsi0, E);
		double x = xa;
		for (int i = 0; i < nStep; i++) {
			RK4Step(x, y, RHS, dx1, nEq);
			x += dx1;

			out << x << "," << y[0] << "," << y[2] << endl;
		}

		// Backward integration
		reset_y(y, nEq, psi1, dpsi1, E);
		x = xb;
		for (int i = 0; i < nStep; i++) {
			RK4Step(x, y, RHS, dx2, nEq);
			x += dx2;

			out << x << "," << y[0] << "," << y[2] << endl;
		}

		E += 1;
	}

	out.close();
}

void reset_y(double y[], const int& nEq, const double& psi, const double& dpsi, const double& E) {
	const double y0[] = {psi, dpsi, E};
	for (int i = 0; i < nEq; i++) {
		y[i] = y0[i];
	}
}

double Res(double E) {
	const double xa = -10.0, xb = 10.0;
	const double xm = 0.5 * (xa + xb) + 1e-1;
	const int nStep = 800;
	const double dx1 = (xm - xa) / nStep;
	const double dx2 = (xm - xb) / nStep;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = -xa * psi0;
	const double psi1 = exp(-0.5 * xb*xb);
	const double dpsi1 = -xb * psi1;
	double y[3];
	const int nEq = 3;

	// Forward integration
	reset_y(y, nEq, psi0, dpsi0, E);
	double x = xa;
	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx1, nEq);
		x += dx1;
	}
	double yL = y[0];
	double dyL = y[1];

	// Backward integration
	reset_y(y, nEq, psi1, dpsi1, E);
	x = xb;
	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx2, nEq);
		x += dx2;
	}
	double yR = y[0];
	double dyR = y[1];

	const double D = 1.0;
	return (dyL * yR - dyR * yL) / D;
}

void ResPlot() {
	std::ofstream out;
	out.open("res_plot.csv");
	if (!out) exit(1);

	double E = 0.0;
	const double dE = 1.0e-3;

	out << "E,Res" << endl;
	while (E < 5.0) {
		out << E << "," << Res(E) << endl;
		E += dE;
	}

	out.close();
}

double ResZero() {
	double roots[8];
	int nRoots;
	find_roots(Res, 0.1, 4.9, 1.0e-8, roots, nRoots);

	for (int i = 0; i < nRoots; i++) {
		cout << "root #" << i << " = " << roots[i] << endl;
	}
	return roots[0];
}

void Final() {
	std::ofstream out;
	out.open("final.csv");
	if (!out) exit(1);

	out << "x,psi,E" << endl;

	const double xa = -10.0, xb = 10.0;
	const double xm = 0.5 * (xa + xb) + 1e-1;
	const int nStep = 800;
	const double dx1 = (xm - xa) / nStep;
	const double dx2 = (xm - xb) / nStep;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = -xa * psi0;
	const double psi1 = exp(-0.5 * xb*xb);
	const double dpsi1 = -xb * psi1;
	const double E = ResZero();
	double y[] = {psi0, dpsi0, E};
	const int nEq = static_cast<int>(sizeof(y) / sizeof(y[0]));

	// Forward integration
	reset_y(y, nEq, psi0, dpsi0, E);
	double x = xa;
	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx1, nEq);
		x += dx1;

		out << x << "," << y[0] << "," << y[2] << endl;
	}

	// Backward integration
	reset_y(y, nEq, psi1, dpsi1, E);
	x = xb;
	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx2, nEq);
		x += dx2;

		out << x << "," << y[0] << "," << y[2] << endl;
	}

	out.close();
}
