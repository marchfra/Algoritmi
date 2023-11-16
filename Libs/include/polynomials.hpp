/**
 * @file polynomials.hpp
 * 
 * @brief Implementation of the polynomial evaluation methods.
 */
#pragma once

#include <vector>

/**
 * @brief Horner method for polynomial evaluation.
 * 
 * Evaluates the polynomial and its derivative using Horner method.
 * 
 * @param[in] x The point at which to evaluate the polynomial.
 * @param[in] a Array of the coefficients. **NOTE**: `a[0]` is the constant term
 *               and `a[degree]` is the coefficient of x^n.
 * @param[in] degree The degree of the polynomial.
 * @param[out] dpdx The value of the derivative.
 * 
 * @returns The value of the polynomial.
 */
double horner_pol(const double&, const double[], const int, double&);

/**
 * @overload
 * 
 * @brief Horner method for polynomial evaluation.
 * 
 * Evaluates the polynomial using Horner method.
 * 
 * @param[in] x The point at which to evaluate the polynomial.
 * @param[in] a Array of the coefficients. **NOTE**: `a[0]` is the constant term
 *               and `a[degree]` is the coefficient of x^n.
 * @param[in] degree The degree of the polynomial.
 * 
 * @returns The value of the polynomial.
 */
double horner_pol(const double&, const double[], const int);

/**
 * @overload
 * 
 * @brief Horner method for polynomial evaluation.
 * 
 * Evaluates the polynomial and its derivative using Horner method.
 * 
 * @param[in] x The point at which to evaluate the polynomial.
 * @param[in] a Vector of the coefficients. **NOTE**: `a[0]` is the constant
 *              term and `a[degree]` is the coefficient of x^n.
 * @param[out] dpdx The value of the derivative.
 * 
 * @returns The value of the polynomial.
 */
double horner_pol(const double&, const std::vector<double>, double&);

/**
 * @overload
 * 
 * @brief Horner method for polynomial evaluation.
 * 
 * Evaluates the polynomial using Horner method.
 * 
 * @param[in] x The point at which to evaluate the polynomial.
 * @param[in] a Vector of the coefficients. **NOTE**: `a[0]` is the constant
 *              term and `a[degree]` is the coefficient of x^n.
 * 
 * @returns The value of the polynomial.
 */
double horner_pol(const double&, const std::vector<double>);
