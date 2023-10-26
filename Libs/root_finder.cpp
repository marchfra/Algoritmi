#include <iostream>
#include <iomanip>

using namespace std;

#include "debug.hpp"
#include "root_finder.hpp"

// =====================================================================================================================
// Bisection method
// =====================================================================================================================

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

			#if DEBUG == TRUE
			cout << setiosflags(ios::scientific);
			cout << "bisection(): k = " << setw(log10(max_ntry) + 1) << k << "; "
				 << "[a, b] = [" << setw(13) << xa << ", " << setw(13) << xb << "]; xm = " << setw(13) << xm << ";\n"
				 << "                   " << setw(log10(max_ntry) + 1) << ""
				 << "err = " << fabs(xb - xa) << "; fm = " << setw(13) << fm << endl;
			cout << resetiosflags(ios::scientific);
			#endif

			// Check convergence
			if (fabs(xb - xa) < xtol || fabs(fm) < ftol || fm == 0.0) {
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

		}

		cout << "! bisection(): too many steps\n" << endl;
		return 1;

	} 

	cout << "! bisection(): initial interval does not contain any root\n" << endl;
	return 2;
}

int bisection(double (*f)(double), double xa, double xb, const double xtol, double& root) {

	// ///////////////////////////////////////////////////////////
	//
	// Overloading of bisection() without ntry and ftol parameters
	//
	// ///////////////////////////////////////////////////////////
	
	int n;
	return bisection(f, xa, xb, xtol, -1.0, root, n);
}

int bisection(double (*f)(double), double xa, double xb, const double xtol, double& root, int& ntry) {

	// /////////////////////////////////////////////////
	//
	// Overloading of bisection() without ftol parameter
	//
	// /////////////////////////////////////////////////
	
	return bisection(f, xa, xb, xtol, -1.0, root, ntry);
}

int bisection(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root) {

	// /////////////////////////////////////////////////
	//
	// Overloading of bisection() without ntry parameter
	//
	// /////////////////////////////////////////////////
	
	int n;
	return bisection(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// False position method
// =====================================================================================================================

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

			#if DEBUG == TRUE
			cout << setiosflags(ios::scientific);
			cout << "false_position(): k = " << setw(log10(max_ntry) + 1) << k << "; "
				 << "[a, b] = [" << setw(13) << xa << ", " << setw(13) << xb << "]; xm = " << setw(13) << xm << "; "
				 << "fm = " << setw(13) << fm << ";\n"
				 << "                        " << setw(log10(max_ntry) + 1) << ""
				 << "err = " << fabs(del) << endl;
			cout << resetiosflags(ios::scientific);
			#endif

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
			if (fabs(del) < xtol || fabs(fm) < ftol || fm == 0.0) {
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

int false_position(double (*f)(double), double xa, double xb, const double xtol, double& root) {

	// ////////////////////////////////////////////////////////////////
	//
	// Overloading of false_position() without ntry and ftol parameters
	//
	// ////////////////////////////////////////////////////////////////
	
	int n;
	return false_position(f, xa, xb, xtol, -1.0, root, n);
}

int false_position(double (*f)(double), double xa, double xb, const double xtol, double& root, int& ntry) {

	// //////////////////////////////////////////////////////
	//
	// Overloading of false_position() without ftol parameter
	//
	// //////////////////////////////////////////////////////
	
	return false_position(f, xa, xb, xtol, -1.0, root, ntry);
}

int false_position(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root) {

	// //////////////////////////////////////////////////////
	//
	// Overloading of false_position() without ntry parameter
	//
	// //////////////////////////////////////////////////////
	
	int n;
	return false_position(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// Secant method
// =====================================================================================================================

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

		#if DEBUG == TRUE
		cout << setiosflags(ios::scientific);
		cout << "secant(): k = " << setw(log10(max_ntry) + 1) << k << "; "
			 << "[a, b] = [" << setw(13) << xa << ", " << setw(13) << xb << "]; err = " << setw(13) << dx << endl;
		cout << resetiosflags(ios::scientific);
		#endif

		// Shift values
		xa = xb;
		fa = fb;
		xb = xb - dx;
		fb = f(xb);

		// Check convergence
		if (fabs(dx) < xtol || fabs(fb) < ftol || fb == 0.0) {
			ntry = k;
			root = xb;
			return 0;
		}
	}

	cout << "! secant(): too many steps\n" << endl;
	return 1;
}

int secant(double (*f)(double), double xa, double xb, const double xtol, double& root) {

	// ////////////////////////////////////////////////////////
	//
	// Overloading of secant() without ntry and ftol parameters
	//
	// ////////////////////////////////////////////////////////
	
	int n;
	return secant(f, xa, xb, xtol, -1.0, root, n);
}

int secant(double (*f)(double), double xa, double xb, const double xtol, double& root, int& ntry) {

	// //////////////////////////////////////////////
	//
	// Overloading of secant() without ftol parameter
	//
	// //////////////////////////////////////////////
	
	return secant(f, xa, xb, xtol, -1.0, root, ntry);
}

int secant(double (*f)(double), double xa, double xb, const double xtol, const double ftol, double& root) {

	// //////////////////////////////////////////////
	//
	// Overloading of secant() without ntry parameter
	//
	// //////////////////////////////////////////////
	
	int n;
	return secant(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// Newton's method
// =====================================================================================================================

int newton(double (*f)(double), double (*dfdx)(double), double xa, double xb, const double xtol, const double ftol, double& root, int& ntry) {

	// /////////////////////////////////////////////////////////////////
	//
	// Find the root of a function f(x) in a given interval [xa, xb]
	// using Newton's method.

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

			#if DEBUG == TRUE
			cout << setiosflags(ios::scientific);
			cout << "newton(): k = " << setw(log10(max_ntry) + 1) << k << "; "
				 << "xc = " << setw(13) << xc << "; dx = " << setw(13) << dx << endl;
			cout << resetiosflags(ios::scientific);
			#endif

			// Check convergence
			if (fabs(dx) < xtol || fabs(fc) < ftol || fc == 0.0) {
				ntry = k;
				root = xc;
				return 0;
			}
		}
	}

	cout << "! newton(): too many steps\n" << endl;
	return 1;
}

int newton(double (*f)(double), double(*dfdx)(double), double xa, double xb, const double xtol, double& root) {

	// ////////////////////////////////////////////////////////
	//
	// Overloading of newton() without ntry and ftol parameters
	//
	// ////////////////////////////////////////////////////////
	
	int n;
	return newton(f, dfdx, xa, xb, xtol, -1.0, root, n);
}

int newton(double (*f)(double), double(*dfdx)(double), double xa, double xb, const double xtol, double& root, int& ntry) {

	// //////////////////////////////////////////////
	//
	// Overloading of newton() without ftol parameter
	//
	// //////////////////////////////////////////////
	
	return newton(f, dfdx, xa, xb, xtol, -1.0, root, ntry);
}

int newton(double (*f)(double), double(*dfdx)(double), double xa, double xb, const double xtol, const double ftol, double& root) {

	// //////////////////////////////////////////////
	//
	// Overloading of newton() without ntry parameter
	//
	// //////////////////////////////////////////////
	
	int n;
	return newton(f, dfdx, xa, xb, xtol, ftol, root, n);
}
