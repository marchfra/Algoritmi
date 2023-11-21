#include "../include/debug.hpp"
#include "../include/root_finder.hpp"

int findRoots(double (*f)(const double& x), double(*dfdx)(const double& x), const double& xa, const double& xb, const double& tol, double roots[], int& nRoots, const int N, const std::string method) {
	if (N > 128) throw std::invalid_argument("N must be less than 128");

	double xL[128], xR[128];

	bracket(f, xa, xb, xL, xR, N, nRoots);

	if (nRoots == 0) {
		std::cerr << "! findRoots(): The interval does not contain any roots:" << std::endl;
		return 2;
	}

	int flag;
	for (int i = 0; i < nRoots; i++) {
		if (method == "bisection") flag = bisection(f, xL[i], xR[i], tol, roots[i]);
		else if (method == "falsePosition") flag = falsePosition(f, xL[i], xR[i], tol, roots[i]);
		else if (method == "secant") flag = secant(f, xL[i], xR[i], tol, roots[i]);
		else if (method == "newton") flag = newton(f, dfdx, xL[i], xR[i], tol, roots[i]);
		else throw std::invalid_argument("Invalid method argument");


		#if DEBUG == TRUE
			std::cout << "roots[" << i << "] = " << roots[i] << std::endl;
		#endif
	}

	std::cout << "Method used: " << method << std::endl;
	return flag;
}

int findRoots(double (*f)(const double& x), const double& xa, const double& xb, const double& tol, double roots[], int& nRoots, const int N, const std::string method) {
	if (method == "newton") throw std::invalid_argument("Newton method isn't available with this prototype.");

	return findRoots(f, nullptr, xa, xb, tol, roots, nRoots, N, method);
}

void bracket(double (*f)(const double& x), const double& xa, const double& xb, double xL[], double xR[], const int& N, int& nRoots) {
	double dx = (xb - xa) / N;
	double xi = xa;
	double xi_plus_one = xi + dx;
	int root_counter = 0;

	double fL = f(xi), fR;
	for (int i = 0; i < N; i++) {
		fR = f(xi_plus_one);
		if (fL == 0.0 || fL * fR < 0) {		// Check if there's a root in [xi, xi_plus_one)
			xL[root_counter] = xi;
			xR[root_counter] = xi_plus_one;
			
			#if DEBUG == TRUE
				std::cout << "found root in [a, b) = [" << xL[root_counter] << ", " << xR[root_counter] << ")" << std::endl;
			#endif
			
			root_counter++;
		}

		// Shift interval
		fL = fR;
		xi = xi_plus_one;
		xi_plus_one += dx;
	}

	nRoots = root_counter;
}

// =====================================================================================================================
// Bisection method
// =====================================================================================================================

int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
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
			std::cout.setf(std::ios::scientific | std::ios::showpos);

			std::cout << "bisection(): k = " << std::setw(log10(max_ntry) + 1) << k << "; err = " << fabs(xb - xa) << std::endl
				<< "                   " << std::setw(log10(max_ntry) + 1) << ""
				<< "xa = " << std::setw(13) << xa << "\t fa = " << fa << "\n"
				<< "                   " << std::setw(log10(max_ntry) + 1) << ""
				<< "xm = " << xm << "\t fm = " << fm << "\n"
				<< "                   " << std::setw(log10(max_ntry) + 1) << ""
				<< "xb = " << xb << "\t fb = " << fb << "\n";

			std::cout << std::resetiosflags(std::ios::scientific | std::ios::showpos);
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

		std::cerr << "! bisection(): too many steps\n" << std::endl;
		ntry = -1;
		root = nan("");
		return 1;

	} 

	std::cerr << "! bisection(): initial interval does not contain any root\n" << std::endl;
	return 2;
}

int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root) {
	int n;
	return bisection(f, xa, xb, xtol, -1.0, root, n);
}

int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
	return bisection(f, xa, xb, xtol, -1.0, root, ntry);
}

int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root) {
	int n;
	return bisection(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// False position method
// =====================================================================================================================

int falsePosition(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
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
			std::cout.setf(std::ios::scientific | std::ios::showpos);
			std::cout << "falsePosition(): k = " << std::setw(log10(max_ntry) + 1) << k << "; "
				 << "[a, b] = [" << xa << ", " << xb << "]; xm = " << xm << "; "
				 << "fm = " << fm << ";\n"
				 << "                        " << std::setw(log10(max_ntry) + 1) << ""
				 << "err = " << fabs(del) << std::endl;
			std::cout << std::resetiosflags(std::ios::scientific | std::ios::showpos);
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

		std::cerr << "! falsePosition(): too many steps\n" << std::endl;
		ntry = -1;
		root = nan("");

		return 1;

	} 

	std::cerr << "! falsePosition(): initial interval does not contain any root\n" << std::endl;
	ntry = -1;
	root = nan("");

	return 2;
}

int falsePosition(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root) {
	int n;
	return falsePosition(f, xa, xb, xtol, -1.0, root, n);
}

int falsePosition(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
	return falsePosition(f, xa, xb, xtol, -1.0, root, ntry);
}

int falsePosition(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root) {
	int n;
	return falsePosition(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// Secant method
// =====================================================================================================================

int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
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
		std::cout.setf(std::ios::scientific | std::ios::showpos);
		std::cout << "secant(): k = " << std::setw(log10(max_ntry) + 1) << k << "; "
			 << "[a, b] = [" << xa << ", " << xb << "]; err = " << dx << std::endl;
		std::cout << std::resetiosflags(std::ios::scientific | std::ios::showpos);
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

	std::cerr << "! secant(): too many steps\n" << std::endl;
	ntry = -1;
	root = nan("");
	return 1;
}

int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root) {
	int n;
	return secant(f, xa, xb, xtol, -1.0, root, n);
}

int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
	return secant(f, xa, xb, xtol, -1.0, root, ntry);
}

int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root) {
	int n;
	return secant(f, xa, xb, xtol, ftol, root, n);
}

// =====================================================================================================================
// Newton's method
// =====================================================================================================================

int newton(double (*f)(const double& x), double (*dfdx)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
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

	// if (fa * fb < 0) {	// Check that interval contains a solution REMOVE FOR EVEN MULTIPLICITY
	for (int k = 1; k <= max_ntry; k++) {
		fc = f(xc);
		dx = fc / dfdx(xc);
		xc -= dx;

		#if DEBUG == TRUE
		std::cout.setf(std::ios::scientific | std::ios::showpos);
		std::cout << "newton(): k = " << std::setw(log10(max_ntry) + 1) << k << "; "
			 << "xc = " << std::setw(13) << xc << "; dx = " << std::setw(13) << dx << std::endl;
		std::cout << std::resetiosflags(std::ios::scientific | std::ios::showpos);
		#endif

		// Check convergence
		if (fabs(dx) < xtol || fabs(fc) < ftol || fc == 0.0) {
			ntry = k;
			root = xc;
			return 0;
		}
	}
	// }

	std::cerr << "! newton(): too many steps\n" << std::endl;
	ntry = -1;
	root = nan("");
	return 1;
}

int newton(double (*f)(const double& x), double(*dfdx)(const double& x), double xa, double xb, const double& xtol, double& root) {
	int n;
	return newton(f, dfdx, xa, xb, xtol, -1.0, root, n);
}

int newton(double (*f)(const double& x), double(*dfdx)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
	return newton(f, dfdx, xa, xb, xtol, -1.0, root, ntry);
}

int newton(double (*f)(const double& x), double(*dfdx)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root) {
	int n;
	return newton(f, dfdx, xa, xb, xtol, ftol, root, n);
}
