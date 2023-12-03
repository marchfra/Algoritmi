#include "../include/polynomials.hpp"

#include "../include/debug.hpp"

double hornerPol(const double& x, const double a[], const int& degree,
                 double& dpdx) {
	double p = a[degree];
	dpdx     = 0.0;
	for (int i = degree - 1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p    = p * x + a[i];
	}
	return p;
}

double hornerPol(const double& x, const double a[], const int& degree) {
	double p = a[degree];
	for (int i = degree - 1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}

double hornerPol(const double& x, const std::vector<double> a, double& dpdx) {
	int degree = a.size();
	double p   = a[degree];
	dpdx       = 0.0;
	for (int i = degree - 1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p    = p * x + a[i];
	}
	return p;
}

double hornerPol(const double& x, const std::vector<double> a) {
	int degree = a.size();
	double p   = a[degree];
	for (int i = degree - 1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}
