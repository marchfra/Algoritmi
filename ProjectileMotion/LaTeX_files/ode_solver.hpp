/**
 * @file    ode_solver.hpp
 *
 * @brief   Implementation of the ODE Solvers step functions.
 *
 * @author  Francesco Marchisotti
 *
 * @date    2024-05-01
 */
#pragma once

#include <iostream>

/**
 * @brief          Runge-Kutta 4 method step.
 *
 * Takes one step in time (or whatever the independent variable is) using
 * Runge-Kutta 4 method. The system of first order ODEs is dY_i/dt = R_i(t, Y).
 *
 * @param[in]      t        The value of time from which to take the step.
 * @param[in, out] Y        Array containing all the dependent variables.
 * @param[in]      RHSFunc  Pointer to the function containing all the Right
 *                          Hand Sides of the system of equations.
 * @param[in]      dt       The step size.
 * @param[in]      neq      The number of equations (ie the number of
 *                          independent variables) in the system.
 */
void rk4Step(const double& t, double Y[],
             void (*RHSFunc)(const double& t, double Y[], double RHS[]),
             const double& dt, const int& neq);
