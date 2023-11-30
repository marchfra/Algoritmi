/**
 * @file quad.hpp
 *
 * @brief      Implementation of the quadrature methods.
 *
 * @author     Francesco Marchisotti
 */
#pragma once

#include <iostream>
#include <cmath>

/**
 * @brief      Rectangular quad method.
 *
 * Integrates a function using the rectangular quadrature method.
 *
 * @param[in]  F     The integrand function.
 * @param[in]  a,b   The lower and upper bound for the integral.
 * @param[in]  n     The number of sub-intervals to use for the integration.
 *
 * @return     The estimate of the integral.
 */
double rectangularQuad(double (*F)(const double& x), const double& a, const double& b, const int& n);

/**
 * @brief      Midpoint quad method.
 *
 * Integrates a function using the midpoint quadrature method.
 *
 * @param[in]  F     The integrand function.
 * @param[in]  a,b   The lower and upper bound for the integral.
 * @param[in]  n     The number of sub-intervals to use for the integration.
 *
 * @return     The estimate of the integral.
 */
double midpointQuad(double (*F)(const double& x), const double& a, const double& b, const int& n);

/**
 * @brief      Trapezoidal quad method.
 *
 * Integrates a function using the trapezoidal quadrature method.
 *
 * @param[in]  F     The integrand function.
 * @param[in]  a,b   The lower and upper bound for the integral.
 * @param[in]  n     The number of sub-intervals to use for the integration.
 *
 * @return     The estimate of the integral.
 */
double trapezoidalQuad(double (*F)(const double& x), const double& a, const double& b, const int& n);

/**
 * @brief      Simpson quad method.
 *
 * Integrates a function using the Simpson quadrature method.
 *
 * @param[in]  F     The integrand function.
 * @param[in]  a,b   The lower and upper bound for the integral.
 * @param[in]  n     The number of sub-intervals to use for the integration. Must
 *                   be even.
 *
 * @return     The estimate of the integral.
 *
 * @throw      std::invalid_argument Thrown if N is even.
 */
double simpsonQuad(double (*F)(const double& x), const double& a, const double& b, const int& n);

/**
 * @brief         Sets the legendre weights and roots.
 *
 * @param[in,out] weights  Array with the legendre weights.
 * @param[in,out] roots    Array with the legendre roots.
 * @param[in]     Ng       Number of gaussian points.
 */
void setLegendreWeightsAndRoots(double weights[], double roots[], const int& Ng);

/**
 * @brief      Gauss-Legendre quad method.
 *
 * Integrates a function using the Gauss-Legendre quadrature method.
 * This method is used with integrals between finite a and b.
 *
 * @param[in]  F     The integrand function.
 * @param[in]  a,b   The lower and upper bound for the integral.
 * @param[in]  N     The number of sub-intervals to use for the integration.
 * @param[in]  Ng    The number of gaussian points to use for each interval.
 *
 * @return     The estimate of the integral.
 */
double gaussLegendreQuad(double (*F)(const double& x), const double& a, const double& b, const int N = 1, const int Ng = 3);

/**
 * @brief      Gauss-Legendre quad method over a rectangle.
 *
 * Integrates a function using the Gauss-Legendre quadrature method
 * over a rectangular domain.
 *
 * @param[in]  F      The integrand function.
 * @param[in]  xa,xb  The lower and upper bound along the x-axis for the integral.
 * @param[in]  ya,yb  The lower and upper bound along the y-axis for the integral.
 * @param[in]  N      The number of sub-intervals to use for the integration.
 * @param[in]  Ng     The number of gaussian points to use for each interval.
 *
 * @return     The estimate of the integral.
 */
double gaussLegendreQuad2D(double (*F)(const double& x, const double& y), const double& xa, const double& xb, const double& ya, const double& yb, const int N = 1, const int Ng = 3);
