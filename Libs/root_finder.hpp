#include <iostream>
#include <iomanip>

int bisection(double (*)(double), double, double, const double, const double, double&, int&);
int bisection(double (*)(double), double, double, const double, double&);
int bisection(double (*)(double), double, double, const double, double&, int&);
int bisection(double (*)(double), double, double, const double, const double, double&);

int false_position(double (*)(double), double, double, const double, const double, double&, int&);
int false_position(double (*)(double), double, double, const double, double&);
int false_position(double (*)(double), double, double, const double, double&, int&);
int false_position(double (*)(double), double, double, const double, const double, double&);

int secant(double (*)(double), double, double, const double, const double, double&, int&);
int secant(double (*)(double), double, double, const double, double&);
int secant(double (*)(double), double, double, const double, double&, int&);
int secant(double (*)(double), double, double, const double, const double, double&);

int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&, int&);
int newton(double (*)(double), double (*)(double), double, double, const double, double&);
int newton(double (*)(double), double (*)(double), double, double, const double, double&, int&);
int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&);
