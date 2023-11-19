/**
 * @file derivative.hpp
 *
 * @brief      Implementation of the numerical differentiation methods.
 * 
 * @author     Francesco Marchisotti
 * 
 * @date       19/11/2023
 */
#pragma once

/**
 * @brief      Forward difference method for first derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double forwardDiff(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Backward difference method for first derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double backwardDiff(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Central difference method for first derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double centralDiff(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Higher order central difference method for first derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double higherDiff(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Forward difference method for second derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double forwardDiff2(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Backward difference method for second derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double backwardDiff2(double (*f)(const double& x), const double& x, const double& h);

/**
 * @brief      Central difference method for second derivative.
 *
 * @param[in]  f     The function to differentiate.
 * @param[in]  x     The value at which to differentiate.
 * @param[in]  h     The step size to use in the differentiation.
 *
 * @return     The value of the derivative.
 */
double centralDiff2(double (*f)(const double& x), const double& x, const double& h);
