/**
 * @file    root_finder.hpp
 *
 * @brief   Implementation of the root finder methods.
 *
 * @author  Paolino Paperino
 *
 * @date    2023-11-19
 */
#pragma once

#include <iomanip>
#include <iostream>

/**
 * @brief      Find the roots of a function in a given interval.
 *
 * Find the roots of a function f(x) in a given interval [xa, xb]
 * using the specified method. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  dfdx    Pointer to the derivative of the function.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[in]  tol     x-tolerance.
 * @param[out] roots   Array with the roots of f(x).
 * @param[out] nRoots  The number of roots found.
 * @param[in]  N       The number of sub-intervals.
 * @param[in]  method  The root finding method. Accepted values are:
 *                     `bisection`, `falsePosition`, `secant`, `newton`.
 *
 * @return     flag
 *
 * @retval     0       Success.
 * @retval     1       Too many steps.
 * @retval     2       Initial interval doesn't contain any root.
 *
 * @throws     std::invalid_argument  Thrown if `N` > 128.
 * @throws     std::invalid_argument  Thrown if `method` is not among the
 *                                    accepted values.
 * @throws     std::runtime_error     Thrown if roots can't be found inside the
 *                                    interval.
 * @throws     std::runtime_error     Thrown if one of the root finders exceeded
 *                                    the maximum number of steps.
 */
int findRoots(double (*f)(const double& x), double (*dfdx)(const double& x),
              const double& xa, const double& xb, const double& tol,
              double roots[], int& nRoots, const int N = 128,
              const std::string method = "newton");

/**
 * @overload
 *
 * @brief      Find the roots of a function in a given interval (Newton's method
 *             not available).
 *
 * Find the roots of a function f(x) in a given interval [xa, xb]
 * using the specified method. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[in]  tol     x-tolerance.
 * @param[out] roots   Array with the roots of f(x).
 * @param[out] nRoots  The number of roots found.
 * @param[in]  N       The number of sub-intervals.
 * @param[in]  method  The root finding method. Accepted values are:
 *                     `bisection`, `falsePosition`, `secant`.
 *
 * @return     flag
 *
 * @retval     0       Success.
 * @retval     1       Too many steps.
 * @retval     2       Initial interval doesn't contain any root.
 *
 * @throws     std::invalid_argument  Thrown if `N` > 128.
 * @throws     std::invalid_argument  Thrown if `method` is not among the
 * accepted values.
 * @throws     std::runtime_error     Thrown if roots can't be found inside the
 *                                    interval.
 * @throws     std::runtime_error     Thrown if one of the root finders exceeded
 *                                    the maximum number of steps.
 */
int findRoots(double (*f)(const double& x), const double& xa, const double& xb,
              const double& tol, double roots[], int& nRoots, const int N = 128,
              const std::string method = "bisection");

/**
 * @brief      Bracket the roots of a function in a given interval [xa, xb].
 *
 * Works by subdividing the interval in a number of sub-intervals
 * and checking if the function changes sign (an odd number of
 * times) over this interval. If it does, then the interval contains
 * (at least) one root.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[out] xL      Array with the lower bound of the sub-interval containing
 *                     a root.
 * @param[out] xR      Array with the upper bound of the sub-interval containing
 *                     a root.
 * @param[in]  N       The number of sub-intervals.
 * @param[out] nRoots  The number of roots found.
 */
void bracket(double (*f)(const double& x), const double& xa, const double& xb,
             double xL[], double xR[], const int& N, int& nRoots);

/**
 * @brief      Find the root of a function f(x) in a given interval [xa, xb]
 *             using secant method.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root  The root of f(x).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x), double xa, double xb,
           const double& xtol, const double& ftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x) in a given interval [xa, xb]
 *             using secant method.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x), double xa, double xb,
           const double& xtol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x) in a given interval [xa, xb]
 *             using secant method.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x), double xa, double xb,
           const double& xtol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x) in a given interval [xa, xb]
 *             using secant method.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root  The root of f(x).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x), double xa, double xb,
           const double& xtol, const double& ftol, double& root);
