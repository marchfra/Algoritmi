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

const static double gT = 2 * M_PI;
const static double gOmega = 2 * M_PI / gT;
const static double gOmega2 = gOmega*gOmega;

void RHS(double Y[], double R[]);
void RHS(const double& t, double Y[], double R[]);

int main() {
	std::ofstream out;
	out.open("../data/data.csv");
	if (!out) exit(5);


	const double h = 0.02 * gT;
	double t = 0.0;
	const double t_max = 10 * gT;
	const int nStep = (t_max - t) / h;
	const double y0[] = {1.0, 0.0};
	double y[2];
	const int neq = static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));

	out << "t,x,v,E,method" << endl;
	for (int i = 0; i < neq; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		rk2Step(t, y, RHS, h, neq);
		t += h;

		double E = 0.5 * (y[1]*y[1] + gOmega2 * y[0]*y[0]);
		out << t << "," << y[0] << "," << y[1] << "," << E << ",RK2" << endl;
	}

	t = 0.0;
	for (int i = 0; i < neq; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		rk4Step(t, y, RHS, h, neq);
		t += h;

		double E = 0.5 * (y[1]*y[1] + gOmega2 * y[0]*y[0]);
		out << t << "," << y[0] << "," << y[1] << "," << E << ",RK4" << endl;
	}

	t = 0.0;
	for (int i = 0; i < neq; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, h, neq);
		t += h;

		double E = 0.5 * (y[1]*y[1] + gOmega2 * y[0]*y[0]);
		out << t << "," << y[0] << "," << y[1] << "," << E << ",pVerlet" << endl;
	}

	t = 0.0;
	for (int i = 0; i < neq; i++) y[i] = y0[i];
	for (int i = 0; i < nStep; i++) {
		vVerlet(t, y, RHS, h, neq);
		t += h;

		double E = 0.5 * (y[1]*y[1] + gOmega2 * y[0]*y[0]);
		out << t << "," << y[0] << "," << y[1] << "," << E << ",vVerlet" << endl;
	}

	out.close();

	return 0;
}

void RHS(double Y[], double R[]) {
	double x = Y[0];
	double v = Y[1];
	R[0] = v;
	R[1] = -gOmega2 * x;
}

void RHS(const double& t, double Y[], double R[]) {
	RHS(Y, R);
}
