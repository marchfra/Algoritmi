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
	double zero;
	double a = -1.0;
	double b =  1.0;
	double tol = 1.0e-7;
	int n = 0.0;

	int flag;

	// for (int i = 0; i < 100; i++) {
	// 	double x = -4.0 + i * (8.0 / 100);
	// 	cout << x << " " << f2(x) << " " << x*x*x - 3*x*x + x + 5 << endl;
	// }

	cout << "--------- Root finders test ---------" << endl;

	flag = bisection(f1, a, b, tol, -1.0, zero, n);
	cout << flag << "; bisection (#" << n << ") =      " << zero << endl;

	flag = false_position(f1, a, b, tol, -1.0, zero, n);
	cout << flag << "; false_position (#" << n << ") = " << zero << endl;

	flag = secant(f1, a, b, tol, -1.0, zero, n);
	cout << flag << "; secant (#" << n << ") =         " << zero << endl;

	flag = newton(f1, df1dx, a, b, tol, -1.0, zero, n);
	cout << flag << "; newton (#" << n << ") =         " << zero << endl;

	cout << "\n------ Polynomial root finding ------" << endl;
	a = -2;
	b =  0;
	tol = 1.0e-8;

	flag = bisection(f2, a, b, tol, -1.0, zero, n);
	cout << flag << "; bisection (#" << n << ") =      " << zero << endl;

	flag = false_position(f2, a, b, tol, -1.0, zero, n);
	cout << flag << "; false_position (#" << n << ") = " << zero << endl;

	flag = secant(f2, a, b, tol, -1.0, zero, n);
	cout << flag << "; secant (#" << n << ") =         " << zero << endl;

	flag = newton(f2, df2dx, a, b, tol, -1.0, zero, n);
	cout << flag << "; newton (#" << n << ") =         " << zero << endl;

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

