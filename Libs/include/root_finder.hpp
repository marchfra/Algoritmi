/**
 * @file root_finder.hpp
 * 
 * @brief Implementation of the root finder methods.
 * 
 * @todo Reformat the functions to not return a flag and instead throw
 *       exceptions.
 */
#pragma once

#include <iostream>
#include <iomanip>

/**
 * @brief Find the roots of a function in a given interval.
 * 
 * Find the roots of a function f(x) in a given interval [xa, xb]
 * using the specified method. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in] f Pointer to the function.
 * @param[in] dfdx Pointer to the derivative of the function.
 * @param[in] xa, xb Interval containing many roots.
 * @param[in] tol x-tolerance.
 * @param[out] roots Array with the roots of f(x).
 * @param[out] n_roots The number of roots found.
 * @param[in] N The number of sub-intervals.
 * @param[in] method The root finding method. Accepted values are: bisection, 
 *            false_position, secant, newton.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int find_roots(double (*)(double), double (*)(double), const double, const double, const double, double[], int&, const int = 128, const std::string = "newton");

/**
 * @overload
 * 
 * @brief Find the roots of a function in a given interval (Newton's method not
 *         available).
 * 
 * Find the roots of a function f(x) in a given interval [xa, xb]
 * using the specified method. Works by first bracketing the roots
 * and then applying the method on every sub-interval.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval containing many roots.
 * @param[in] tol x-tolerance.
 * @param[out] roots Array with the roots of f(x).
 * @param[out] n_roots The number of roots found.
 * @param[in] N The number of sub-intervals.
 * @param[in] method The root finding method. Accepted values are: bisection, 
 *            false_position, secant.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 * 
 * @throws std::invalid_argument Thrown if `method` is `"newton"`.
 */
int find_roots(double (*)(double), const double, const double, const double, double[], int&, const int = 128, const std::string = "bisection");


/**
 * @brief Bracket the roots of a function in a given interval [xa, xb].
 * 
 * Works by subdividing the interval in a number of sub-intervals and checking 
 * if the function changes sign (an odd number of times) over this interval.
 * If it does, then the interval contains (at least) one root.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval containing many roots.
 * @param[out] xL,xR Arrays with the left and right endpoints of the interval
 *                   containing a root.
 * @param[in] N The number of sub-intervals.
 * @param[out] n_roots The number of roots found.
 */
void bracket(double (*)(double), const double, const double, double[], double[], const int, int&);


/**
 * @brief Find the root of a funciont in a given interval using bisection 
 *        method.
 * 
 * Find the root of a function f(x) in a given interval [xa, xb]
 * using bisection method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int bisection(double (*)(double), double, double, const double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using bisection 
 *        method.
 * 
 * @param[in] f Pointer to the function
 * @param[in] xa, xb Interval (containing the root)
 * @param[in] xtol x-tolerance
 * @param[out] root The root of f(x)
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int bisection(double (*)(double), double, double, const double, double&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using bisection 
 *        method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int bisection(double (*)(double), double, double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using bisection 
 *        method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int bisection(double (*)(double), double, double, const double, const double, double&);

/**
 * @brief Find the root of a funciont in a given interval using false position 
 *        method.
 * 
 * Find the root of a function f(x) in a given interval [xa, xb]
 * using false position method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int false_position(double (*)(double), double, double, const double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using false position 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int false_position(double (*)(double), double, double, const double, double&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using false position 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int false_position(double (*)(double), double, double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using false position 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 * @retval 2 Initial interval doesn't contain any root.
 */
int false_position(double (*)(double), double, double, const double, const double, double&);


/**
 * @brief Find the root of a funciont in a given interval using the secant 
 *        method.
 * 
 * Find the root of a function f(x) in a given interval [xa, xb]
 * using the secant method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int secant(double (*)(double), double, double, const double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using the secant 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int secant(double (*)(double), double, double, const double, double&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using the secant 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int secant(double (*)(double), double, double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using the secant 
 *        method.
 *
 * @param[in] f Pointer to the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int secant(double (*)(double), double, double, const double, const double, double&);

/**
 * @brief Find the root of a funciont in a given interval using Newton-Raphson 
 *        method.
 * 
 * Find the root of a function f(x) in a given interval [xa, xb]
 * using Newton-Rhapson method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] dfdx Pointer to the derivative of the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using Newton-Raphson 
 *        method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] dfdx Pointer to the derivative of the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int newton(double (*)(double), double (*)(double), double, double, const double, double&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using Newton-Raphson 
 *        method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] dfdx Pointer to the derivative of the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[out] root The root of f(x).
 * @param[out] ntry The number of iterations achieved.
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int newton(double (*)(double), double (*)(double), double, double, const double, double&, int&);

/**
 * @overload
 * 
 * @brief Find the root of a funciont in a given interval using Newton-Raphson 
 *        method.
 * 
 * @param[in] f Pointer to the function.
 * @param[in] dfdx Pointer to the derivative of the function.
 * @param[in] xa, xb Interval (containing the root).
 * @param[in] xtol x-tolerance.
 * @param[in] ftol f(x)-tolerance: the values of f(x) that are considered 0.
 * @param[out] root The root of f(x).
 *
 * @retval 0 Success.
 * @retval 1 Too many steps.
 */
int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&);
