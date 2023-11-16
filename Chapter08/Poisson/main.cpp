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

double rho(double r);

void RHS(double r, double Y[], double R[]);

double Residual(double s);

double Exact(double r);

int main() {
	double s = 0.0;
	double r = 0.0;
	const double r_max = 20;
	const int nStep = 1000;
	const double dr = (r_max - r) / nStep;

	// Shooting
	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);

	out << "r,phi,phi',Phi,s" << endl;
	while (s <= 1.0) {
		double y[] = {0.0, s};
		const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

		r = 0.0;
		for (int i = 0; i < nStep; i++) {
			RK4Step(r, y, RHS, dr, nEq);
			r += dr;

			out << r << "," << y[0] << "," << y[1] << "," << y[0]/r << "," << s << endl;
		}

		s += 0.2;
	}

	out.close();

	// Residual plot
	out.open("residual.csv");
	if (!out) exit(1);

	out << "s,Res" << endl;
	s = 0.0;
	const double s_max = 5.0;
	const int nS = 10;
	const double ds = (s_max - s) / nS;

	for (int i = 0; i < nS; i++) {
		out << s << "," << Residual(s) << endl;
		s += ds;
	}

	out.close();

	// Find the root of Residual(s)
	out.open("final.csv");
	if (!out) exit(1);

	double xa = 0.0;
	double xb = 1.0;
	double tol = 1.0e-8;
	double s_root;
	bisection(Residual, xa, xb, tol, s_root);

	cout << "s_root = " << s_root << endl;

	double y_final[] = {0.0, s_root};
	const int nEq = static_cast<int>(sizeof(y_final)) / static_cast<int>(sizeof(y_final[0]));


	out << "r,phi,Phi,exact,Exact,s" << endl;
	r = 0.0;
	for (int i = 0; i < nStep; i++) {
		RK4Step(r, y_final, RHS, dr, nEq);
		r += dr;

		out << r << "," << y_final[0] << "," << y_final[0]/r << "," << Exact(r) << "," << Exact(r)/r << "," << s_root << endl;
	}


	out.close();
	return 0;
}

double rho(double r) {
	return exp(-r) / (8.0 * M_PI);
}

void RHS(double r, double Y[], double R[]) {
	R[0] = Y[1];
	R[1] = -4.0 * M_PI * r * rho(r);
}

double Residual(double s) {
	double r = 0.0;
	const double r_max = 20;
	const int nStep = 1000;
	const double dr = (r_max - r) / nStep;

	double y[] = {0.0, s};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	for (int i = 0; i < nStep; i++) {
		RK4Step(r, y, RHS, dr, nEq);
		r += dr;
	}

	return y[0] - 1.0;		// phi(b, s) - 1
}

double Exact(double r) {
	return 1.0 - 0.5 * (r + 2) * exp(-r);
}
