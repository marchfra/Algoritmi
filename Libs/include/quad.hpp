/**
 * @file quad.hpp
 * 
 * @brief Implementation of the quadrature methods.
 * 
 * @todo Reformat the functions to throw exceptions.
 */
#pragma once

#include <iostream>
#include <cmath>

/**
 * @brief Rectangular quad method.
 * 
 * Integrates a function using the rectangular quadrature method.
 * 
 * @param[in] F The integrand function.
 * @param[in] a, b The lower and upper bound for the integral.
 * @param[in] n The number of sub-intervals to use for the integration.
 * 
 * @returns The estimate of the integral.
 */
double rectangular_quad(double (*)(double), const double, const double, const int);

/**
 * @brief Midpoint quad method.
 * 
 * Integrates a function using the midpoint quadrature method.
 * 
 * @param[in] F The integrand function.
 * @param[in] a, b The lower and upper bound for the integral.
 * @param[in] n The number of sub-intervals to use for the integration.
 * 
 * @returns The estimate of the integral.
 */
double midpoint_quad(double (*)(double), const double, const double, const int);

/**
 * @brief Trapezoidal quad method.
 * 
 * Integrates a function using the trapezoidal quadrature method.
 * 
 * @param[in] F The integrand function.
 * @param[in] a, b The lower and upper bound for the integral.
 * @param[in] n The number of sub-intervals to use for the integration.
 * 
 * @returns The estimate of the integral.
 */
double trapezoidal_quad(double (*)(double), const double, const double, const int);

/**
 * @brief Simpson quad method.
 * 
 * Integrates a function using the Simpson quadrature method.
 * 
 * @param[in] F The integrand function.
 * @param[in] a, b The lower and upper bound for the integral.
 * @param[in] n The number of sub-intervals to use for the integration. Must be 
 *              even.
 * 
 * @returns The estimate of the integral.
 */
double simpson_quad(double (*)(double), const double, const double, const int);

/**
 * @brief Gauss-Legendre quad method.
 * 
 * Integrates a function using the Gauss-Legendre quadrature method.
 * This method is used with integrals between finite a and b.
 * 
 * @param[in] F The integrand function.
 * @param[in] a, b The lower and upper bound for the integral.
 * @param[in] N The number of sub-intervals to use for the integration.
 * @param[in] Ng The number of gaussian points to use for each interval.
 * 
 * @returns The estimate of the integral.
 */
double gauss_legendre_quad(double (*)(double), const double, const double, const int = 1, const int = 3);

/**
 * @brief Gauss-Legendre quad method over a rectangle.
 * 
 * Integrates a function using the Gauss-Legendre quadrature method over a 
 * rectangular domain.
 * 
 * @param[in] F The integrand function.
 * @param[in] xa, xb The lower and upper bound along the x-axis for the integral.
 * @param[in] ya, yb The lower and upper bound along the y-axis for the integral.
 * @param[in] N The number of sub-intervals to use for the integration.
 * @param[in] Ng The number of gaussian points to use for each interval.
 * 
 * @returns The estimate of the integral.
 */
double gauss_legendre_quad2D(double (*)(double, double), const double, const double, const double, const double, const int = 1, const int = 3);
