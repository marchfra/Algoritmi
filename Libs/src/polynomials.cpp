#include <iostream>

#include "../include/debug.hpp"
#include "../include/polynomials.hpp"

double horner_pol(const double& x, const double a[], const int degree, double& dpdx) {

	// //////////////////////////////////////////////////////////////
	//
	// Evaluates the polynomial Σ_n a_n x^n and its derivative in x
	// 
	// x        [in]  : point at which to evaluate the polynomial
	// a[]      [in]  : array of the coefficients
	//                  NOTE: a[0] is the constant term and a[degree]
	//                        is the coefficient of x^n
	// degree   [in]  : the degree of the polynomial
	// dpdx     [out] : the value of the derivative at x
	//
	// Returns: the value of the polynomial at x
	//
	// Last modified: 26 Oct 2023
	//
	// //////////////////////////////////////////////////////////////

	double p = a[degree];
	dpdx = 0.0;
	for (int i = degree-1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p = p * x + a[i];
	}
	return p;
}

double horner_pol(const double& x, const double a[], const int degree) {

	// /////////////////////////////////////////////////////////////
	//
	// Overloading of horner_pol() without derivative implementation
	//
	// /////////////////////////////////////////////////////////////

	double p = a[degree];
	for (int i = degree-1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}

double horner_pol(const double& x, const double a[]) {

	// /////////////////////////////////////////////////////////////
	//
	// Overloading of horner_pol() without derivative implementation
	// and with automatic degree detection
	//
	// /////////////////////////////////////////////////////////////

	std::cout << "! horner_pol(): this overloading doesn't work" << std::endl;
	exit(1);

	const int degree = static_cast<int>(sizeof(*a)) / static_cast<int>(sizeof(a[0]));

	double p = a[degree];
	for (int i = degree-1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}
