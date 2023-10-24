#include "polynomials.hpp"

double horner_pol(const double x, const double a[], const int degree) {
	double p = a[degree];
	for (int i = degree-1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}

double horner_pol(const double x, const double a[], const int degree, double& dpdx) {
	double p = a[degree];
	dpdx = 0.0;
	for (int i = degree-1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p = p * x + a[i];
	}
	return p;
}
