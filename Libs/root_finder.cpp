#include <iostream>
#include <iomanip>

using namespace std;

#include "root_finder.hpp"

// #define DEBUG FALSE

int bisection(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root, int& ntry) {

	// //////////////////////////////////////////////////////////////
	//
	// Find the root of a function f(x) in a given interval [xa, xb]
	// using bisection method.

	// *f       [in]  : pointer to the function
	// xa, xb   [in]  : initial interval (containing the root)
	// xtol     [in]  : x-tolerance
	// ftol     [in]  : f(x)-tolerance
	// root     [out] : the root of f(x)
	// ntry     [out] : number of iterations achieved
	//
	// On output the function returns 0 (success), 1 (too many steps)
	// or 2 (initial interval doesn't contain any root).
	//
	// Last modified: 26 Oct 2023
	//
	// //////////////////////////////////////////////////////////////


	int max_ntry = 128;
	double fa = f(xa);
	double fb = f(xb);
	double xm, fm;

	// Handle fa, fb = 0
	if (fa == 0.0) {
		ntry = 0;
		root = xa;
		return 0;
	} else if (fb == 0.0) {
		ntry = 0;
		root = xb;
		return 0;
	}


	if (fa * fb < 0) {		// Necessary condition
		for (int k = 1; k <= max_ntry; k++) {
			// Function midpoint
			xm = 0.5 * (xa + xb);
			fm = f(xm);

			// Check convergence
			if (fabs(xb - xa) < xtol || fabs(fm) < ftol) {
				ntry = k;
				root = xm;
				return 0;
			}

			// Redefine interval
			if (fm * fa < 0) {
				xb = xm;
				fb = fm;
			} else {
				xa = xm;
				fa = fm;
			}

			#if DEBUG == TRUE
			cout << scientific << "bisection(): k = " << k << "; [a, b] = [" << xa << ", " << xb << "]; xm = " << xm << "; err = " << fabs(xb - xa) << "; fm = " << fm << endl;
			#endif

		}

		cout << "! bisection(): too many steps\n" << endl;
		return 1;

	} 

	cout << "! bisection(): initial interval does not contain any root\n" << endl;
	return 2;
}

int false_position(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root, int& ntry) {

	// //////////////////////////////////////////////////////////////
	//
	// Find the root of a function f(x) in a given interval [xa, xb]
	// using false_position method.

	// *f       [in]  : pointer to the function
	// xa, xb   [in]  : initial interval (containing the root)
	// xtol     [in]  : x-tolerance
	// ftol     [in]  : f(x)-tolerance
	// root     [out] : the root of f(x)
	// ntry     [out] : number of iterations achieved
	//
	// On output the function returns 0 (success), 1 (too many steps)
	// or 2 (initial interval doesn't contain any root).
	//
	// Last modified: 26 Oct 2023
	//
	// //////////////////////////////////////////////////////////////


	int max_ntry = 128;
	double fa = f(xa);
	double fb = f(xb);
	double xm, fm;
	double del = xb - xa;

	// Handle fa, fb = 0
	if (fa == 0.0) {
		ntry = 0;
		root = xa;
		return 0;
	} else if (fb == 0.0) {
		ntry = 0;
		root = xb;
		return 0;
	}


	if (fa * fb < 0) {		// Necessary condition
		for (int k = 1; k <= max_ntry; k++) {
			// Linear intersection
			xm = (xa * fb - xb * fa) / (fb - fa);
			fm = f(xm);

			// Redefine interval
			if (fm * fa < 0) {
				del = xb - xm;
				xb = xm;
				fb = fm;
			} else {
				del = xm - xa;
				xa = xm;
				fa = fm;
			}

			// Check convergence
			if (fabs(del) < xtol || fabs(fm) < ftol) {
				ntry = k;
				root = xm;
				return 0;
			}
		}

		cout << "! false_position(): too many steps\n" << endl;
		return 1;

	} 

	cout << "! false_position(): initial interval does not contain any root\n" << endl;
	return 2;
}

int secant(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root, int& ntry) {

	// /////////////////////////////////////////////////////////////////
	//
	// Find the root of a function f(x) in a given interval [xa, xb]
	// using secant method.

	// *f       [in]    pointer to the function
	// xa, xb   [in]    initial interval (containing the root)
	// xtol     [in]    x-tolerance
	// ftol     [in]    f(x)-tolerance
	// root     [out]   the root of f(x)
	// ntry     [out]   number of iterations achieved
	//
	// On output the function returns 0 (success) or 1 (too many steps).
	//
	// Last modified: 26 Oct 2023
	//
	// /////////////////////////////////////////////////////////////////

	int max_ntry = 64;
	double fa = f(xa);
	double fb = f(xb);
	double dx = xb - xa;

	// Handle fa, fb = 0
	if (fa == 0.0) {
		ntry = 0;
		root = xa;
		return 0;
	} else if (fb == 0.0) {
		ntry = 0;
		root = xb;
		return 0;
	}

	for (int k = 1; k <= max_ntry; k++) {
		dx = fb * (xb - xa) / (fb - fa);	// Compute increment

		#if DEBUG
		cout << scientific << "secant(): k = " << k << "; xa = " << xa << "; xb = " << xb << "; dx = " << dx << endl;
		#endif

		// Shift values
		xa = xb;
		fa = fb;
		xb = xb - dx;
		fb = f(xb);

		// Check convergence
		if (fabs(dx) < xtol || fabs(fb) < ftol) {
			ntry = k;
			root = xb;
			return 0;
		}
	}

	cout << "! secant(): too many steps\n" << endl;
	return 1;
}

int newton(double (*f)(double), double (*dfdx)(double), double xa, double xb, const double xtol, const double ftol, double& root, int& ntry) {

	// /////////////////////////////////////////////////////////////////
	//
	// Find the root of a function f(x) in a given interval [xa, xb]
	// using secant method.

	// *f       [in]  : pointer to the function
	// *dfdx    [in]  : pointer to the derivative of the function
	// xa, xb   [in]  : initial interval (containing the root)
	// xtol     [in]  : x-tolerance
	// ftol     [in]  : f(x)-tolerance
	// root     [out] : the root of f(x)
	// ntry     [out] : number of iterations achieved
	//
	// On output the function returns 0 (success) or 1 (too many steps).
	//
	// Last modified: 26 Oct 2023
	//
	// /////////////////////////////////////////////////////////////////

	int max_ntry = 16;
	double fa = f(xa);
	double fb = f(xb);
	double fc, dx;
	double xc = 0.5 * (xa + xb);

	// Handle fa, fb = 0
	if (fa == 0.0) {
		ntry = 0;
		root = xa;
		return 0;
	} else if (fb == 0.0) {
		ntry = 0;
		root = xb;
		return 0;
	}

	if (fa * fb < 0) {	// Check that interval contains a solution REMOVE FOR EVEN MULTIPLICITY
		for (int k = 1; k <= max_ntry; k++) {
			fc = f(xc);
			dx = fc / dfdx(xc);
			xc -= dx;

			#if DEBUG
			cout << scientific << "newton(): k = " << k << "; xc = " << xc << "; dx = " << dx << endl;
			#endif

			// Check convergence
			if (fabs(dx) < xtol || fabs(fc) < ftol) {
				ntry = k;
				root = xc;
				return 0;
			}
		}
	}

	cout << "! newton(): too many steps\n" << endl;
	return 1;
}
