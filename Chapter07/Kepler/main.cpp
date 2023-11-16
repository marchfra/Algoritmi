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

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

using std::setw;

void R(double t, double Y[], double R[]);

void ellipse_test();

void N_orbits();

int main() {
	// ellipse_test();
	N_orbits();

	return 0;
}

void R(double t, double Y[], double R[]) {
	double x = Y[0];
	double y = Y[1];
	double vx = Y[2];
	double vy = Y[3];
	const double r = sqrt(x*x + y*y);
	const double r3 = r*r*r;

	R[0] = vx;
	R[1] = vy;
	R[2] = -x / r3;
	R[3] = -y / r3;
}

void ellipse_test() {
	std::ofstream out;
	out.open("ellipse_test.csv");
	if (!out) exit(1);

	const double x0 = 4.0, y0 = 0.0;
	const double r0 = sqrt(x0*x0 + y0*y0);
	const double alpha = 0.4;
	const double vx0 = 0.0, vy0 = sqrt(alpha / r0);
	double y[] = {x0, y0, vx0, vy0};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	const double t0 = 0.0;
	const double tend = 100.0;
	double t = t0;
	const int nStep = 1000;
	const double dt = (tend - t0) / nStep;

	// cout << "  t         x(t)          y(t)           vx(t)          vy(t)" << endl;
	out << "t,x,y,vx,vy" << endl;
	for (int i = 0; i < nStep; i++) {
		RK4Step(t, y, R, dt, nEq);
		t += dt;
		// cout << std::resetiosflags(std::ios::scientific);
		// cout << setw(4) << t << "  ";
		// cout.setf(std::ios::scientific);
		// cout << setw(13) << y[0] << "  " << setw(13) << y[1] << "  ";
		// cout << setw(13) << y[2] << "  " << setw(13) << y[3] << endl;
		out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << endl;
	}

	out.close();
}

void N_orbits() {
	std::ofstream out;
	out.open("N_orbits.csv");
	if (!out) exit(1);

	const double x0 = 4.0, y0 = 0.0;
	const double r0 = sqrt(x0*x0 + y0*y0);
	const double alpha = 0.3;
	const double vx0 = 0.0, vy0 = sqrt(alpha / r0);
	double y[] = {x0, y0, vx0, vy0};
	const int nEq = static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	const double t0 = 0.0;
	double t = t0;
	const int nMaxStep = 10000;
	const double dTheta = 0.1;
	const int nMaxOrbits = 10;
	int nOrbit = 0;
	// double prev_vx = y[2];
	double prev_r2 = r0*r0;
	double pprev_r2 = prev_r2;
	const double E0 = 0.5 * (vx0*vx0 + vy0*vy0) - 1 / r0;
	double E = E0;
	const double L0 = y0 * vx0 - x0 * vy0;
	double L = L0;

	out << "t,x,y,vx,vy,E,L" << endl;
	out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << E << "," << L << endl;
	double dt;
	for (int i = 0; i < nMaxStep; i++) {
		double r2 = y[0]*y[0] + y[1]*y[1];	// Compute r^2
		double v2 = y[2]*y[2] + y[3]*y[3];	// Compute v^2

		// Adaptive dt: maintain dTheta ~constant
		// dt = dTheta * r / v_t, v_t = tangential velocity (not easy)
		// I use v insead of v_t, thus underestimating dt. Oh no. Anyway...
		dt = dTheta * sqrt(r2 / v2);

		RK4Step(t, y, R, dt, nEq);
		t += dt;

		E = 0.5 * v2 - 1 / sqrt(r2);
		L = y[1] * y[2] - y[0] * y[3];
		out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << E << "," << L << endl;
		// if (prev_vx * y[2] < 0.0) {
		if (r2 < prev_r2 && prev_r2 > pprev_r2) {
			nOrbit++;
			if (2*nOrbit == nMaxOrbits) {
				cout << "Reached " << nMaxOrbits << " orbits" << endl;
				cout << "Made " << i << " steps" << endl;
				break;
			}
		}
		// prev_vx = y[2];
		pprev_r2 = prev_r2;
		prev_r2 = r2;
	}

	out.close();
}
