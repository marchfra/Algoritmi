#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/quad.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/polynomials.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/root_finder.hpp"

using namespace std;

void test_methods(double (*)(double), double (*)(double), const double, const double, const double, double&, int&);

double f1(double);
double df1dx(double);

double f2(double);
double df2dx(double);

double f3(double);
double df3dx(double);

int main() {
	double root;
	double tol = 1.0e-7;
	int n;

	cout << "--------- Root finders test ---------" << endl;
	double a = -1.0;
	double b =  1.0;

	test_methods(f1, df1dx, a, b, tol, root, n);

	cout << "\n------ Polynomial root finding ------" << endl;
	tol = 1.0e-8;
	a = -5;
	b =  0;

	cout << "[a, b] = [" << a << ", " << b << "]" << endl;

	test_methods(f2, df2dx, a, b, tol, root, n);

	cout << endl;
	a = -2;
	b =  0;
	cout << "[a, b] = [" << a << ", " << b << "]" << endl;

	test_methods(f2, df2dx, a, b, tol, root, n);

	cout << "\n------ Ugly thing root finding ------" << endl;
	tol = 1.0e-7;
	a = 0;
	b = 2;

	test_methods(f3, df3dx, a, b, tol, root, n);

	return 0;
}

void test_methods(double (*f)(double), double (*dfdx)(double), const double a, const double b, const double tol, double& root, int& n) {
	bisection(f, a, b, tol, root, n);
	cout << "bisection (#" << n << ") = " << root << endl;

	false_position(f, a, b, tol, root, n);
	cout << "false_pos (#" << n << ") = " << root << endl;

	secant(f, a, b, tol, root, n);
	cout << "secant (#" << n << ") = " << root << endl;

	newton(f, dfdx, a, b, tol, root, n);
	cout << "newton (#" << n << ") = " << root << endl;
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
}

double df2dx(double x) {
	double a[] = {5.0, 1.0, -3.0, 1.0};
	double dpdx;
	horner_pol(x, a, 4, dpdx);
	return dpdx;
}

double f3(double x) {
	return exp(1 / (x + 0.5)) - (3 + 2 * x) / (1 + x);
}

double df3dx(double x) {
	return 1 / ((1 + x)*(1 + x)) - exp(1 / (x + 0.5)) / ((x + 0.5)*(x + 0.5));
}

