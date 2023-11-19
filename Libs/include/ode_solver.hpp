/**
 * @file ode_solver.hpp
 *
 * @brief      Implementation of the ODE Solvers step functions.
 * 
 * @todo       Fix RK4Step.\n Improve pVerlet using pointers to improve
 *             readability.\n Fix vVerlet and rewrite using pointers.
 * @author     Francesco Marchisotti
 * 
 * @date       19/11/2023
 */
#pragma once

#include <iostream>

/**
 * @brief          Euler method step.
 *
 * Takes one step in time (or whatever the independent variable is) using Euler
 * method. The system of first order ODEs is dY_i/dt = R_i(t, Y).
 *
 * @param[in]      t        The value of time from which to take the step.
 * @param[in, out] Y        Array containing all the dependent variables.
 * @param[in]      RHSFunc  Pointer to the function containing all the Right
 *                          Hand Sides of the system of equations.
 * @param[in]      dt       The step size.
 * @param[in]      neq      The number of equations (ie the number of
 *                          independent variables) in the system.
 */
void eulerStep(const double& t, double Y[], void (*RHSFunc)(const double& t, double Y[], double RHS[]), const double& dt, const int& neq);

/**
 * @brief          Runge-Kutta 2 method step.
 *
 * Takes one step in time (or whatever the independent variable is) using
 * Runge-Kutta 2 method. The system of first order ODEs is dY_i/dt = R_i(t, Y).
 *
 * @param[in]      t        The value of time from which to take the step.
 * @param[in, out] Y        Array containing all the dependent variables.
 * @param[in]      RHSFunc  Pointer to the function containing all the Right
 *                          Hand Sides of the system of equations.
 * @param[in]      dt       The step size.
 * @param[in]      neq      The number of equations (ie the number of
 *                          independent variables) in the system.
 * @param[in]      rule     If `false` use Midpoint method (default). If `true`
 *                          use Modified Euler method.
 */
void rk2Step(const double& t, double Y[], void (*RHSFunc)(const double& t, double Y[], double RHS[]), const double& dt, const int& neq, bool rule = false);

/**
 * @brief          Runge-Kutta 4 method step.
 *
 * Takes one step in time (or whatever the independent variable is) using 
 * Runge-Kutta 4 method. The system of first order ODEs is dY_i/dt = R_i(t, Y).
 *
 * @bug            This method converges rather than diverging.
 *
 * @param[in]      t        The value of time from which to take the step.
 * @param[in, out] Y        Array containing all the dependent variables.
 * @param[in]      RHSFunc  Pointer to the function containing all the Right
 *                          Hand Sides of the system of equations.
 * @param[in]      dt       The step size.
 * @param[in]      neq      The number of equations (ie the number of
 *                          independent variables) in the system.
 */
void rk4Step(const double& t, double Y[], void (*RHSFunc)(const double& t, double Y[], double RHS[]), const double& dt, const int& neq);

/**
 * @brief          Position Verlet method step.
 *
 * Takes one step in time using Position Verlet method. The system of first order 
 * ODEs is dY_i/dt = R_i(Y), where all equations are derived from second order
 * ODEs in the form F = ma.
 *
 * @param[in]      t        The value of time from which to take the step.
 * @param[in, out] Y        Array containing all the dependent variables.
 * @param[in]      RHSFunc  Pointer to the function containing all the Right
 *                          Hand Sides of the system of equations.
 * @param[in]      dt       The step size.
 * @param[in]      neq      The number of equations (ie the number of
 *                          independent variables) in the system. Must be even.
 */
void pVerlet(const double& t, double Y[], void (*RHSFunc)(double Y[], double RHS[]), const double& dt, const int& neq);

/**
 * @brief      Velocity Verlet method step.
 *
 * Takes one step in time using Velocity Verlet method. The system of first order
 * ODEs is dY_i/dt = R_i(Y), where all equations are derived from second order
 * ODEs in the form F = ma.
 *
 * @param[in]  t        The value of time from which to take the step.
 * @param[in]  Y        Array containing all the dependent variables.
 * @param[in]  RHSFunc  Pointer to the function containing all the Right Hand
 *                      Sides of the system of equations.
 * @param[in]  dt       The step size.
 * @param[in]  neq      The number of equations (ie the number of independent
 *                      variables) in the system. Must be even.
 */
void vVerlet(const double& t, double Y[], void (*RHSFunc)(double Y[], double RHS[]), const double& dt, const int& neq);
