/**
 * @file polynomials.hpp
 *
 * @brief      Implementation of the polynomial evaluation methods.
 * 
 * @author     Francesco Marchisotti.
 * 
 * @date       19/11/2023
 */
#pragma once

#include <vector>

/**
 * @brief      Horner method for polynomial evaluation.
 *
 * Evaluates the polynomial and its derivative using Horner method.
 *
 * @param[in]  x       The point at which to evaluate the polynomial.
 * @param[in]  a       Array of the coefficients. **NOTE**: `a[0]` is the
 *                     constant term and `a[degree]` is the coefficient of x^n.
 * @param[in]  degree  The degree of the polynomial.
 * @param[out] dpdx    The value of the derivative.
 *
 * @return     The value of the polynomial.
 */
double hornerPol(const double& x, const double a[], const int& degree, double& dpdx);

/**
 * @overload
 *
 * @brief      Horner method for polynomial evaluation.
 *
 * Evaluates the polynomial using Horner method.
 *
 * @param[in]  x       The point at which to evaluate the polynomial.
 * @param[in]  a       Array of the coefficients. **NOTE**: `a[0]` is the
 *                     constant term and `a[degree]` is the coefficient of x^n.
 * @param[in]  degree  The degree of the polynomial.
 *
 * @return     The value of the polynomial.
 */
double hornerPol(const double& x, const double a[], const int& degree);

/**
 * @overload
 *
 * @brief      Horner method for polynomial evaluation.
 *
 * Evaluates the polynomial and its derivative using Horner method.
 *
 * @param[in]  x     The point at which to evaluate the polynomial.
 * @param[in]  a     Vector of the coefficients. **NOTE**: `a[0]` is the
 *                   constant term and `a[degree]` is the coefficient of x^n.
 * @param[out] dpdx  The value of the derivative.
 *
 * @return     The value of the polynomial.
 */
double hornerPol(const double& x, const std::vector<double> a, double& dpdx);

/**
 * @overload
 *
 * @brief      Horner method for polynomial evaluation.
 *
 * Evaluates the polynomial using Horner method.
 *
 * @param[in]  x     The point at which to evaluate the polynomial.
 * @param[in]  a     Vector of the coefficients. **NOTE**: `a[0]` is the
 *                   constant term and `a[degree]` is the coefficient of x^n.
 *
 * @return     The value of the polynomial.
 */
double hornerPol(const double& x, const std::vector<double> a);
