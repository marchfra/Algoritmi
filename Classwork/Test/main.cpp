#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

double f(double x);

double pol(double x);
double dpdx(double x);

int main() {
	cout << "Hello World!" << endl;

	cout << "\n---------- Test swap ----------" << endl;
	int nA = 0, nB = 5;
	float fA = 3.14f, fB = 9.81f;
	double dA = 14.7, dB = 1.62;
	std::string strA = "test1", strB = "pie";

	cout << "nA = " << nA << "; nB = " << nB << endl;
	cout << "Swap(nA, nB)" << endl;
	swap(nA, nB);
	cout << "nA = " << nA << "; nB = " << nB << endl << endl;

	cout << "fA = " << fA << "; fB = " << fB << endl;
	cout << "Swap(fA, fB)" << endl;
	swap(fA, fB);
	cout << "fA = " << fA << "; fB = " << fB << endl << endl;

	cout << "dA = " << dA << "; dB = " << dB << endl;
	cout << "Swap(dA, dB)" << endl;
	swap(dA, dB);
	cout << "dA = " << dA << "; dB = " << dB << endl << endl;

	cout << "strA = " << strA << "; strB = " << strB << endl;
	cout << "Swap(strA, strB)" << endl;
	swap(strA, strB);
	cout << "strA = " << strA << "; strB = " << strB << endl << endl;

	cout << "\n---------- Test quad ----------" << endl;
	double dXA = 0.0, dXB = 5.0;
	int nN = 6;
	cout << "Integral from 0 to 5 of f(x) = x^2 - 2" << endl;
	cout << "Rectangle: " << rectangularQuad(f, dXA, dXB, nN) << endl;
	cout << "Trapezoid: " << trapezoidal_quad(f, dXA, dXB, nN) << endl;
	cout << "Simpson:   " << simpson_quad(f, dXA, dXB, nN) << endl;
	cout << "Gauss:     " << gauss_legendre_quad(f, dXA, dXB, nN) << endl;

	cout << "\n---------- Test poly ----------" << endl;
	cout << "Integral from 0 to 5 of f(x) = x^2 - 2" << endl;
	cout << "Rectangle: " << rectangularQuad(pol, dXA, dXB, nN) << endl;
	cout << "Trapezoid: " << trapezoidal_quad(pol, dXA, dXB, nN) << endl;
	cout << "Simpson:   " << simpson_quad(pol, dXA, dXB, nN) << endl;
	cout << "Gauss:     " << gauss_legendre_quad(pol, dXA, dXB, nN) << endl;

	cout << "\n---------- Test root ----------" << endl;
	double roots[15];
	int n_roots;
	const int n_intervals = 10;
	const double tol = 1.0e-10;
	
	find_roots(pol, -10.0, 10.0, tol, roots, n_roots, n_intervals);
	for (int i = 0; i < n_roots; i++) {
		cout << "root[" << i << "] = " << roots[i] << endl;
	}

	find_roots(pol, -10.0, 10.0, tol, roots, n_roots, n_intervals, "false_position");
	for (int i = 0; i < n_roots; i++) {
		cout << "root[" << i << "] = " << roots[i] << endl;
	}

	find_roots(pol, -10.0, 10.0, tol, roots, n_roots, n_intervals, "secant");
	for (int i = 0; i < n_roots; i++) {
		cout << "root[" << i << "] = " << roots[i] << endl;
	}

	find_roots(pol, dpdx, -10.0, 10.0, tol, roots, n_roots, n_intervals);
	for (int i = 0; i < n_roots; i++) {
		cout << "root[" << i << "] = " << roots[i] << endl;
	}

	return 0;
}

double f(double x) {
	return x*x - 2.0;
}

double pol(double x) {
	double a[] = {-2.0, 0.0, 1.0};
	return horner_pol(x, a, 2);
}

double dpdx(double x) {
	double a[] = {-2.0, 0.0, 1.0};
	double dp = 0.0;
	horner_pol(x, a, 2, dp);
	return dp;
}
