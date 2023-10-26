#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/quad.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/polynomials.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/root_finder.hpp"

using namespace std;

double f1(double);
double df1dx(double);

double f2(double);
double df2dx(double);

int main() {
	double root;
	double tol = 1.0e-7;
	int n;

	cout << "--------- Root finders test ---------" << endl;
	double a = -1.0;
	double b =  1.0;

	bisection(f1, a, b, tol, root, n);
	cout << "bisection (#" << n << ") = " << root << endl;

	false_position(f1, a, b, tol, root, n);
	cout << "false_position (#" << n << ") = " << root << endl;

	secant(f1, a, b, tol, root, n);
	cout << "secant (#" << n << ") = " << root << endl;

	newton(f1, df1dx, a, b, tol, root, n);
	cout << "newton (#" << n << ") = " << root << endl;

	cout << "\n------ Polynomial root finding ------" << endl;
	tol = 1.0e-8;
	a = -5;
	b =  0;

	cout << "Interval [a, b] = [" << a << ", " << b << "]" << endl;

	bisection(f2, a, b, tol, root, n);
	cout << "bisection (#" << n << ") = " << root << endl;

	false_position(f2, a, b, tol, root, n);
	cout << "false_position (#" << n << ") = " << root << endl;

	secant(f2, a, b, tol, root, n);
	cout << "secant (#" << n << ") = " << root << endl;

	newton(f2, df2dx, a, b, tol, root, n);
	cout << "newton (#" << n << ") = " << root << endl;

	cout << endl;
	a = -2;
	b =  0;
	cout << "Interval [a, b] = [" << a << ", " << b << "]" << endl;

	bisection(f2, a, b, tol, root, n);
	cout << "bisection (#" << n << ") = " << root << endl;

	false_position(f2, a, b, tol, root, n);
	cout << "false_position (#" << n << ") = " << root << endl;

	secant(f2, a, b, tol, root, n);
	cout << "secant (#" << n << ") = " << root << endl;

	newton(f2, df2dx, a, b, tol, root, n);
	cout << "newton (#" << n << ") = " << root << endl;

	return 0;
}

double f1(double x) {
	return exp(-x) - x;
}

double df1dx(double x) {
	return -exp(-x) - 1;
}

double f2(double x) {
	double a[] = {5.0, 1.0, -3.0, 1.0};
	return horner_pol(x, a, 4);
	// return x*x*x - 3*x*x + x + 5;
}

double df2dx(double x) {
	double a[] = {5.0, 1.0, -3.0, 1.0};
	double dpdx;
	horner_pol(x, a, 4, dpdx);
	return dpdx;
	// return 3*x*x - 6*x + 1;
}

