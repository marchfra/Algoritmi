#include "../include/derivative.hpp"

#include "../include/debug.hpp"

double forwardDiff(double (*f)(const double &x), const double &x,
                   const double &h) {
  double fxph = f(x + h);
  double fx = f(x);

  return (fxph - fx) / h;
}

double backwardDiff(double (*f)(const double &x), const double &x,
                    const double &h) {
  double fx = f(x);
  double fxmh = f(x - h);

  return (fx - fxmh) / h;
}

double centralDiff(double (*f)(const double &x), const double &x,
                   const double &h) {
  double fxph = f(x + h);
  double fxmh = f(x - h);

  return (fxph - fxmh) / (2.0 * h);
}

double higherDiff(double (*f)(const double &x), const double &x,
                  const double &h) {
  double fxpph = f(x + 2.0 * h);
  double fxph = f(x + h);
  double fxmh = f(x - h);
  double fxmmh = f(x - 2.0 * h);

  return (fxmmh - 8.0 * fxmh + 8.0 * fxph - fxpph) / (12.0 * h);
}

double forwardDiff2(double (*f)(const double &x), const double &x,
                    const double &h) {
  double fxpph = f(x + 2.0 * h);
  double fxph = f(x + h);
  double fx = f(x);

  return (fxpph - 2.0 * fxph + fx) / (h * h);
}

double backwardDiff2(double (*f)(const double &x), const double &x,
                     const double &h) {
  double fx = f(x);
  double fxmh = f(x - h);
  double fxmmh = f(x - 2.0 * h);

  return (fxmmh - 2.0 * fxmh + fx) / (h * h);
}

double centralDiff2(double (*f)(const double &x), const double &x,
                    const double &h) {
  double fxph = f(x + h);
  double fx = f(x);
  double fxmh = f(x - h);

  return (fxph - 2.0 * fx + fxmh) / (h * h);
}
