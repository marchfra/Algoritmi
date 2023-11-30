#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <climits>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;

void R(const double& t, double Y[], double R[]);

void convergence_test();

int main() {
	const double t0 = 0.0;
	const double te = 2'000.0 * M_PI;
	double y0[] = {1.0, 0.0};
	double y[2];
	const int nStep = 20'000;
	const double h = te / nStep;

	// cout << "      t              x(t)            y(t)" << endl; //         abs_err        rel_err" << endl;
	// cout.setf(std::ios::scientific | std::ios::showpos);
	// double t = t0;
	// while (t <= te) {
	// 	EulerStep(t, y, R, h, 2);
	// 	t += h;
	// 	cout << t << "   " << y[0] << "   " << y[1] << endl; // << abs_err << "   " << rel_err << endl;
	// }


	std::ofstream out;
	out.open("data/data.csv");
	if (!out) exit(1);

	out << "t,x,y,method" << endl;
	out.setf(std::ios::scientific | std::ios::showpos);
	// Euler integration
	double t = t0;
	for (int i = 0; i < 2; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		eulerStep(t, y, R, h, 2);
		t += h;
		out << t << "," << y[0] << "," << y[1] << ",Euler" << endl;
	}

	// RK2 integration
	t = t0;
	for (int i = 0; i < 2; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		rk2Step(t, y, R, h, 2);
		t += h;
		out << t << "," << y[0] << "," << y[1] << ",RK2" << endl;
	}

	// RK4 integration
	t = t0;
	for (int i = 0; i < 2; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		rk4Step(t, y, R, h, 2);
		t += h;
		out << t << "," << y[0] << "," << y[1] << ",RK4" << endl;
	}

	out.close();

	convergence_test();
	return 0;
}

void R(const double& t, double Y[], double R[]) {
	double x = Y[0];
	double y = Y[1];
	R[0] =  y;
	R[1] = -x;
}

void convergence_test() {
	const double t0 = 0.0;
	const double te = 3.0;
	double y0[] = {1.0, 0.0};
	double y[2];
	int nStep = 4;
	double h = te / nStep;

	std::ofstream out;
	out.open("data/convergence.csv");

	out << "dt,err,method" << endl;
	out.setf(std::ios::scientific | std::ios::showpos);
	while (nStep <= 2048) {
		h = te / nStep;

		// Euler integration
		double t = t0;
		double err;
		for (int i = 0; i < 2; i++) y[i] = y0[i];
		for (int i = 0; i < nStep; i++) {
			eulerStep(t, y, R, h, 2);
			t += h;
		}
		err = fabs(y[0] - cos(te));
		out << h << "," << err << ",Euler" << endl;

		// RK2 integration
		t = t0;
		for (int i = 0; i < 2; i++) y[i] = y0[i];
		for (int i = 0; i < nStep; i++) {
			rk2Step(t, y, R, h, 2);
			t += h;
		}
		err = fabs(y[0] - cos(te));
		out << h << "," << err << ",RK2" << endl;

		// RK4 integration
		t = t0;
		for (int i = 0; i < 2; i++) y[i] = y0[i];
		for (int i = 0; i < nStep; i++) {
			rk4Step(t, y, R, h, 2);
			t += h;
		}
		err = fabs(y[0] - cos(te));
		out << h << "," << err << ",RK4" << endl;

		nStep *= 2;
	}

	out.close();
}
