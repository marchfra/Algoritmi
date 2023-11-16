/**
 * @file derivative.hpp
 * 
 * @brief Implementation of the numerical differentiation methods.
 * 
 * @todo Reformat the code to be consistent with other libraries in naming
 *       scheme.
 */
#pragma once

/**
 * @brief Forward difference method for first derivative.
 * 
 * Forward difference method for first derivative (first order accurate).
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double ForwardDiff(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Backward difference method for first derivative.
 * 
 * Backward difference method for first derivative (first order accurate).
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double BackwardDiff(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Central difference method for first derivative.
 * 
 * Central difference method for first derivative (second order accurate).
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double CentralDiff(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Higher order central difference method for first derivative.
 * 
 * Central difference method for first derivative (fourth order accurate).
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double HigherDiff(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Forward difference method for second derivative.
 * 
 * Forward difference method for second derivative.
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double ForwardDiff2(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Backward difference method for second derivative.
 * 
 * Backward difference method for second derivative.
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double BackwardDiff2(const double (*f)(const double), const double dX, const double dH);

/**
 * @brief Central difference method for second derivative.
 * 
 * Central difference method for second derivative.
 * 
 * @param[in] f The function to differentiate.
 * @param[in] dX The value at which to differentiate.
 * @param[in] dH The step size to use in the differentiation.
 * 
 * @returns The value of the derivative.
 */
double CentralDiff2(const double (*f)(const double), const double dX, const double dH);
