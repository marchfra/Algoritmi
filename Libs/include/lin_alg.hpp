/**
 * @file lin_alg.hpp
 *
 * @brief      This file implements linear algorithms.
 *
 * @author     Francesco Marchisotti
 *
 * @todo       Write tests
 *
 * @date       24/11/2023
 */

#pragma once

/**
 * @brief      Prints a matrix.
 *
 * @param[in]  M      `n x n` square matrix.
 * @param[in]  n      Size of M.
 * @param[in]  width  Minimum width of the printing area.
 */
void printMatrix(const double **M, const int& n, const int width=8);

/**
 * @brief      Prints a vector.
 *
 * @param[in]  v      The vector.
 * @param[in]  nRows  Size of the vector.
 */
void printVector(const double v[], const int& nRows);

/**
 * @brief         Swap two rows of a linear system in matrix form.
 *
 * @param[in,out] M      The `n x n` coefficient matrix.
 * @param[in,out] v      The constant vector.
 * @param[in]     nRows  The number of equations in the system.
 * @param[in]     i      The first row to be swapped.
 * @param[in]     j      The second row to be swapped.
 */
void swapRowsLinSystem(double **M, double v[], const int& nRows, const int& i, const int& j);

/**
 * @brief         Gauss reduction algorithm.
 *
 * @param[in,out] M      The coefficient matrix.
 * @param[in,out] v      The constant matrix.
 * @param[in]     nRows  The number of equations.
 */
void gaussElim(double **M, double v[], const int& nRows);

/**
 * @brief         Partial pivoting algorithm.
 *
 * @param[in,out] M            The coefficient matrix.
 * @param[in,out] v            The constant vector.
 * @param[in]     nRows        The number of rows of M and v.
 * @param[in]     current_row  The current row.
 */
void partialPivoting(double **M, double v[], const int& nRows, const int& current_row);

/**
 * @brief         Solve a linear system of equations in matrix form.
 *
 * @param[in,out] M     The coefficient matrix.
 * @param[in]     v     The constant vector.
 * @param[out]    x     The variable vector.
 * @param[in]     nEqs  The number of equations.
 */
void solveLinSystem(double **M, double v[], double x[], const int& nEqs);
