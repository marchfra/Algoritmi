#pragma once

#include <iostream>
#include <cmath>

double rectangular_quad(double (*)(double), const double, const double, const int);

double midpoint_quad(double (*)(double), const double, const double, const int);

double trapezoidal_quad(double (*)(double), const double, const double, const int);

double simpson_quad(double (*)(double), const double, const double, const int);

double gauss_legendre_quad(double (*)(double), const double, const double, const int = 1, const int = 3);

// Integrate over a rectangle
double gauss_legendre_quad2D(double (*)(double, double), const double, const double, const double, const double, const int = 1, const int = 3);
