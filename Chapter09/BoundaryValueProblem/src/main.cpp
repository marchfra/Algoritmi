#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

double rFunc(const double& x);

int main() {
	const int nPoints = 32;

	// Boundary conditions
	const double xL = 0.0;
	const double xR = 1.0;
	const double yL = 1.0;
	const double yR = 0.9;

	double y[nPoints];
	const double dx = (xR - xL) / (nPoints - 1);

	linearBVP(rFunc, y, xL, xR, yL, yR, nPoints);

	std::ofstream out;
	out.open("../data/data.csv");
	if (!out) exit(5);
	out << "x,y" << endl;
	for (int i = 0; i < nPoints; i++) {
		out << xL + i * dx << "," << y[i] << endl;
	}
	out.close();
	return 0;
}

double rFunc(const double& x) {
	return 1.0;
}
