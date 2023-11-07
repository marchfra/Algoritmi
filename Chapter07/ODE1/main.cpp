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

void R(double t, double Y[], double R[]);

double exact(double t);

int main() {
	const double t0 = 0.0;
	const double te = 3.0;
	double h  = 1.0;
	double y[] = {1.0};

	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);

	cout.setf(std::ios::scientific);
	double abs_err, rel_err;
	// cout << "     t              y(t)         abs_err        rel_err" << endl;
	// cout << t << "   " << y[0] << "   " << abs_err << "   " << rel_err << endl;
	out << "t,y,abs_err,rel_err,h" << endl;
	while (h >= 0.001) {
		abs_err = rel_err = 0.0;
		double t = t0;
		y[0] = 1.0;
		out << t << "," << y[0] << "," << abs_err << "," << rel_err << "," << h << endl;
		while (t <= te) {
			EulerStep(t, y, R, h, 1);
			double ex = exact(t);
			abs_err = fabs(y[0] - ex);
			rel_err = abs_err / ex;
			// cout << t << "   " << y[0] << "   " << abs_err << "   " << rel_err << endl;
			t += h;
			out << t << "," << y[0] << "," << abs_err << "," << rel_err << "," << h << endl;
		}
		h /= 2.0;
	}

	out.close();

	return 0;
}

void R(double t, double Y[], double R[]) {
	double y = Y[0];
	R[0] = -t * y;
}

double exact(double t) {
	return exp(-0.5 * t*t);
}
