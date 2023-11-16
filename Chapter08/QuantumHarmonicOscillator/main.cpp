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

int main() {
	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);

	out << "x,psi" << endl;

	const double xa = -10.0, xb = 10.0;
	const int nStep = 800;
	const double dx = (xb - xa) / nStep;
	double x = xa;

	const double psi0 = exp(-0.5 * xa*xa);
	const double dpsi0 = xa * psi0;
	const double E0 = 0.5;
	double y[] = {psi0, dpsi0, E0};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	for (int i = 0; i < nStep; i++) {
		RK4Step(x, y, RHS, dx, nEq);
		x += dx;

		out << x << "," << y[0] << endl;
	}

	out.close();
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
