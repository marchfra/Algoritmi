#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

const static double b = 1;

void RHS(double Y[], double R[]);

int main() {
	std::ofstream out;
	out.open("../data/data.csv");
	if (!out) exit(5);

	// Integration range
	double t = 0.0;
	const double t_max = 10.0;
	const double dt = 0.01;
	const int nStep = (t_max - t) / dt;

	// Initialisation
	const double theta[] = {
		40.0 / 180 * M_PI,
		45.0 / 180 * M_PI,
		50.0 / 180 * M_PI};
	const double v0 = 1.0;
	const double y0[] = {0.0, 0.0, v0 * cos(theta[0]), v0 * sin(theta[0])};
	const double y1[] = {0.0, 0.0, v0 * cos(theta[1]), v0 * sin(theta[1])};
	const double y2[] = {0.0, 0.0, v0 * cos(theta[2]), v0 * sin(theta[2])};
	double y[4];
	const int nEq = static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));

	out << "t,x,y,u,v,theta,drag_coeff" << endl;

	for (int i = 0; i < nEq; i++) y[i] = y0[i];
	out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[0] << "," << b << endl;
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[0] << endl;

		if (y[1] < 0.0) break;
	}

	t = 0.0;
	for (int i = 0; i < nEq; i++) y[i] = y1[i];
	out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[1] << "," << b << endl;
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[1] << "," << b << endl;

		if (y[1] < 0.0) break;
	}

	t = 0.0;
	for (int i = 0; i < nEq; i++) y[i] = y2[i];
	out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[2] << "," << b << endl;
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta[2] << endl;

		if (y[1] < 0.0) break;
	}


	out.close();
	return 0;
}

void RHS(double Y[], double R[]) {
	// double x = Y[0];
	// double y = Y[1];
	double u = Y[2];
	double v = Y[3];

	double mod_v = sqrt(u*u + v*v);

	R[0] = u;
	R[1] = v;
	R[2] = -b * u * mod_v;
	R[3] = -b * 1.0 - u * mod_v;
}
