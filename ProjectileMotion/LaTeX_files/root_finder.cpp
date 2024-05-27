#include "../include/root_finder.hpp"

int findRoots(double (*f)(const double &x), double (*dfdx)(const double &x),
              const double &xa, const double &xb, const double &tol,
              double roots[], int &nRoots, const int N,
              const std::string method) {
  if (N > 128) throw std::invalid_argument("N must be less than 128");

  double xL[128], xR[128];

  bracket(f, xa, xb, xL, xR, N, nRoots);

  if (nRoots == 0) {
    throw std::runtime_error(
      "The supplied interval does not contain any roots.");
  }

  for (int i = 0; i < nRoots; i++)
    if (method == "bisection") bisection(f, xL[i], xR[i], tol, roots[i]);
    else if (method == "falsePosition")
      falsePosition(f, xL[i], xR[i], tol, roots[i]);
    else if (method == "secant") secant(f, xL[i], xR[i], tol, roots[i]);
    else if (method == "newton") newton(f, dfdx, xL[i], xR[i], tol, roots[i]);
    else throw std::invalid_argument("Invalid method argument.");

  return 0;
}

int findRoots(double (*f)(const double &x), const double &xa, const double &xb,
              const double &tol, double roots[], int &nRoots, const int N,
              const std::string method) {
  if (method == "newton")
    throw std::invalid_argument(
      "Newton method isn't available with this prototype.");

  return findRoots(f, nullptr, xa, xb, tol, roots, nRoots, N, method);
}

void bracket(double (*f)(const double &x), const double &xa, const double &xb,
             double xL[], double xR[], const int &N, int &nRoots) {
  double dx          = (xb - xa) / N;
  double xi          = xa;
  double xi_plus_one = xi + dx;
  int root_counter   = 0;

  double fL = f(xi), fR;
  for (int i = 0; i < N; i++) {
    fR = f(xi_plus_one);
    if (fL == 0.0 ||
        fL * fR < 0) {  // Check if there's a root in [xi, xi_plus_one)
      xL[root_counter] = xi;
      xR[root_counter] = xi_plus_one;

      root_counter++;
    }

    // Shift interval
    fL = fR;
    xi = xi_plus_one;
    xi_plus_one += dx;
  }

  nRoots = root_counter;
}

// =============================================================================
// Secant method
// =============================================================================

int secant(double (*f)(const double &x), double xa, double xb,
           const double &xtol, const double &ftol, double &root, int &ntry) {
  int max_ntry = 64;
  double fa    = f(xa);
  double fb    = f(xb);
  double dx    = xb - xa;

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
    dx = fb * (xb - xa) / (fb - fa);  // Compute increment

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

  ntry = -1;
  root = nan("");
  throw std::runtime_error("Maximum number of steps exceeded.");
  return 1;
}

int secant(double (*f)(const double &x), double xa, double xb,
           const double &xtol, double &root) {
  int n;
  return secant(f, xa, xb, xtol, -1.0, root, n);
}

int secant(double (*f)(const double &x), double xa, double xb,
           const double &xtol, double &root, int &ntry) {
  return secant(f, xa, xb, xtol, -1.0, root, ntry);
}

int secant(double (*f)(const double &x), double xa, double xb,
           const double &xtol, const double &ftol, double &root) {
  int n;
  return secant(f, xa, xb, xtol, ftol, root, n);
}
