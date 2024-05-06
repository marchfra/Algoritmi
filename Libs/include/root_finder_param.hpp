/**
 * @file       root_finder_param.hpp
 *
 * @brief      Implementation of the root finder methods.
 *
 * @author     Francesco Marchisotti
 *
 * @date       06/05/2024
 */
#pragma once

#include <iostream>
#include <iomanip>

/**
 * @overload
 *
 * @brief      Find the roots of a function in a given interval.
 *
 * Find the roots of a function f(x, k) in a given interval [xa, xb]
 * using the specified method. k is a parameter of the function and is simply passed to the function. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  dfdx    Pointer to the derivative of the function.
 * @param[in]  k       Additional function parameter.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[in]  tol     x-tolerance.
 * @param[out] roots   Array with the roots of f(x, k).
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
int findRoots(double (*f)(const double& x, const double& k), double (*dfdx)(const double& x, const double& k), const double& k, const double& xa, const double& xb, const double& tol, double roots[], int& nRoots, const int N = 128, const std::string method = "newton");

/**
 * @overload
 *
 * @brief      Find the roots of a function in a given interval (Newton's method
 *             not available).
 *
 * Find the roots of a function f(x, k) in a given interval [xa, xb]
 * using the specified method. k is a parameter of the function and is simply passed to the function. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  k       Additional function parameter.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[in]  tol     x-tolerance.
 * @param[out] roots   Array with the roots of f(x, k).
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
 * @throws     std::invalid_argument  Thrown if `method` is not among the accepted
 *                                    values.
 * @throws     std::runtime_error     Thrown if roots can't be found inside the
 *                                    interval.
 * @throws     std::runtime_error     Thrown if one of the root finders exceeded
 *                                    the maximum number of steps.
 */
int findRoots(double (*f)(const double& x, const double& k), const double& k, const double& xa, const double& xb, const double& tol, double roots[], int& nRoots, const int N = 128, const std::string method = "bisection");

/**
 * @brief      Bracket the roots of a function in a given interval [xa, xb].
 *
 * Works by subdividing the interval in a number of sub-intervals
 * and checking if the function changes sign (an odd number of
 * times) over this interval. If it does, then the interval contains
 * (at least) one root.
 *
 * @param[in]  f       Pointer to the function.
 * @param[in]  k       Additional function parameter.
 * @param[in]  xa      Lower bound of the interval.
 * @param[in]  xb      Upper bound of the interval.
 * @param[out] xL      Array with the lower bound of the sub-interval containing
 *                     a root.
 * @param[out] xR      Array with the upper bound of the sub-interval containing
 *                     a root.
 * @param[in]  N       The number of sub-intervals.
 * @param[out] nRoots  The number of roots found.
 */
void bracket(double (*f)(const double& x, const double& k), const double& k, const double& xa, const double& xb, double xL[], double xR[], const int& N, int& nRoots);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using bisection method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int bisection(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using bisection method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int bisection(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using bisection method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int bisection(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using bisection method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int bisection(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using false position method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int falsePosition(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using false position method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int falsePosition(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using false position method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 * @param[out] ntry  The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int falsePosition(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using false position method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 * @retval     2     Initial interval doesn't contain any root.
 *
 * @throws     std::runtime_error  Thrown if roots can't be found inside the
 *                                 interval.
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int falsePosition(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root);

/**
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using secant method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
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
int secant(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using secant method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using secant method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
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
int secant(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using secant method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[in]  ftol  f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[out] root  The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 */
int secant(double (*f)(const double& x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, double& root);

/**
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using Newton-Raphson method. k is simply passed to the function.
 *
 * @param[in]  f      Pointer to the function.
 * @param[in]  dfdx   Pointer to the derivative of the function.
 * @param[in]  k      Additional function parameter.
 * @param[in]  xa     Lower bound of the interval.
 * @param[in]  xb     Upper bound of the interval.
 * @param[in]  xtol   x-tolerance.
 * @param[in]  ftol   f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[in]  dftol  f'(x, k)-tolerance: protects against obvious divergence. Pass
 *                    -1 to disable.
 * @param[out] root   The root of f(x, k).
 * @param[out] ntry   The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0      Success.
 * @retval     1      Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 * @throws     std::runtime_error  Thrown if the derivative becomes too small
 *                                 (divergence protection).
 */
int newton(double (*f)(const double& x, const double& k), double (*dfdx)(const double &x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, const double& dftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using Newton-Raphson method. k is simply passed to the function.
 *
 * @param[in]  f     Pointer to the function.
 * @param[in]  dfdx  Pointer to the derivative of the function.
 * @param[in]  k     Additional function parameter.
 * @param[in]  xa    Lower bound of the interval.
 * @param[in]  xb    Upper bound of the interval.
 * @param[in]  xtol  x-tolerance.
 * @param[out] root  The root of f(x, k).
 *
 * @note       Tolerance on f'(x, k) is 1e-3. To disable tolerance run `newton()`
 *             [?/4] with `dftol = -1`.
 *
 * @return     flag
 *
 * @retval     0     Success.
 * @retval     1     Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 * @throws     std::runtime_error  Thrown if the derivative becomes too small
 *                                 (divergence protection).
 */
int newton(double (*f)(const double& x, const double& k), double (*dfdx)(const double &x, const double& k), const double& k, double xa, double xb, const double& xtol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using Newton-Raphson method. k is simply passed to the function.
 *
 * @param[in]  f      Pointer to the function.
 * @param[in]  dfdx   Pointer to the derivative of the function.
 * @param[in]  k      Additional function parameter.
 * @param[in]  xa     Lower bound of the interval.
 * @param[in]  xb     Upper bound of the interval.
 * @param[in]  xtol   x-tolerance.
 * @param[in]  dftol  f'(x, k)-tolerance: protects against obvious divergence. Pass
 *                    -1 to disable.
 * @param[out] root   The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0      Success.
 * @retval     1      Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 * @throws     std::runtime_error  Thrown if the derivative becomes too small
 *                                 (divergence protection).
 */
int newton(double (*f)(const double& x, const double& k), double (*dfdx)(const double &x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& dftol, double& root);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using Newton-Raphson method. k is simply passed to the function.
 *
 * @param[in]  f      Pointer to the function.
 * @param[in]  dfdx   Pointer to the derivative of the function.
 * @param[in]  k      Additional function parameter.
 * @param[in]  xa     Lower bound of the interval.
 * @param[in]  xb     Upper bound of the interval.
 * @param[in]  xtol   x-tolerance.
 * @param[in]  dftol  f'(x, k)-tolerance: protects against obvious divergence. Pass
 *                    -1 to disable.
 * @param[out] root   The root of f(x, k).
 * @param[out] ntry   The number of iterations achieved.
 *
 * @return     flag
 *
 * @retval     0      Success.
 * @retval     1      Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 * @throws     std::runtime_error  Thrown if the derivative becomes too small
 *                                 (divergence protection).
 */
int newton(double (*f)(const double& x, const double& k), double (*dfdx)(const double &x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& dftol, double& root, int& ntry);

/**
 * @overload
 *
 * @brief      Find the root of a function f(x, k) in a given interval [xa, xb]
 *             using Newton-Raphson method. k is simply passed to the function.
 *
 * @param[in]  f      Pointer to the function.
 * @param[in]  dfdx   Pointer to the derivative of the function.
 * @param[in]  k      Additional function parameter.
 * @param[in]  xa     Lower bound of the interval.
 * @param[in]  xb     Upper bound of the interval.
 * @param[in]  xtol   x-tolerance.
 * @param[in]  ftol   f(x, k)-tolerance: the values of f(x, k) that are considered 0.
 * @param[in]  dftol  f'(x, k)-tolerance: protects against obvious divergence. Pass
 *                    -1 to disable.
 * @param[out] root   The root of f(x, k).
 *
 * @return     flag
 *
 * @retval     0      Success.
 * @retval     1      Too many steps.
 *
 * @throws     std::runtime_error  Thrown if the maximum number of steps is
 *                                 exceeded.
 * @throws     std::runtime_error  Thrown if the derivative becomes too small
 *                                 (divergence protection).
 */
int newton(double (*f)(const double& x, const double& k), double (*dfdx)(const double &x, const double& k), const double& k, double xa, double xb, const double& xtol, const double& ftol, const double& dftol, double& root);
