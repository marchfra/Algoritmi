/**
 * @file lin_alg.hpp
 *
 * @brief      This file implements linear algorithms.
 *
 * @author     Francesco Marchisotti
 *
 * @date       24/11/2023
 */

#pragma once

#include <iostream>
#include <iomanip>

#include "../include/swap.hpp"

/**
 * @brief      Prints a matrix.
 *
 * @param[in]  M      `n x n` square matrix.
 * @param[in]  n      Size of M.
 * @param[in]  width  Minimum width of the printing area.
 *
 * @tparam     T      Type of the elements of the matrix.
 */
template <class T>
void printMatrix(T **M, const int& n, const int width = 8) {
	for (int i = 0; i < n; i++) {
		std::cout << std::string(width + 1, '-');
	}
	std::cout << std::endl;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << std::setw(width) << M[i][j] << " ";
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < n; i++) {
		std::cout << std::string(width + 1, '-');
	}
	std::cout << std::endl;
}

/**
 * @brief      Prints a vector.
 *
 * @param[in]  v      The vector.
 * @param[in]  nRows  Size of the vector.
 *
 * @tparam     T      Type of the elements of the vector.
 */
template <class T>
void printVector(T v[], const int& nRows) {
	std::cout << "{";
	for (int i = 0; i < nRows; i++) {
		std::cout << v[i];
		if (i != nRows - 1) std::cout << ", ";
	}
	std::cout << "}" << std::endl;
}

/**
 * @brief         Swap two rows of a linear system in matrix form.
 *
 * @param[in,out] M      The `n x n` coefficient matrix.
 * @param[in,out] v      The constant vector.
 * @param[in]     nRows  The number of equations in the system.
 * @param[in]     i      The first row to be swapped.
 * @param[in]     j      The second row to be swapped.
 *
 * @tparam        T      Type of the elements in `M` and `v`.
 */
template <class T>
void swapRowsLinSystem(T **M, T v[], const int& nRows, const int& i, const int& j) {
	for (int k = 0; k < nRows; k++) {
		swap(M[i][k], M[j][k]);
	}

	swap(v[i], v[j]);
}

/**
 * @brief         Partial pivoting algorithm.
 *
 * @param[in,out] M            The coefficient matrix.
 * @param[in,out] v            The constant vector.
 * @param[in]     nRows        The number of rows of M and v.
 * @param[in]     current_row  The current row.
 *
 * @tparam        T            Type of the elements in `M` and `v`.
 */
template <class T>
void partialPivoting(T **M, T v[], const int& nRows, const int& current_row) {
	T max = fabs(M[current_row][current_row]), maxI = current_row;
	for (int i = current_row; i < nRows; i++) {	// Loop over rows under current row
		T g = fabs(M[i][current_row]);
		if (g > max) {
			max = g;
			maxI = i;
		}
	}
	if (maxI != current_row) {
		swapRowsLinSystem(M, v, nRows, maxI, current_row);
	}
}

/**
 * @brief         Gauss reduction algorithm.
 *
 * @param[in,out] M      The coefficient matrix.
 * @param[in,out] v      The constant matrix.
 * @param[in]     nRows  The number of equations.
 *
 * @tparam        T      Type of the elements in `M` and `v`.
 */
template <class T>
void gaussElim(T **M, T v[], const int& nRows) {
	// Reduce M to upper triangular form
	for (int k = 0; k < nRows - 1; k++) {
		partialPivoting(M, v, nRows, k);
		for (int i = k + 1; i < nRows; i++) {
			T g = M[i][k] / M[k][k];
			for (int j = k + 1; j < nRows; j++) M[i][j] -= g * M[k][j];
			M[i][k] = 0.0;
			v[i] -= g * v[k];
		}
	}
}

/**
 * @brief         Solve a linear system of equations in matrix form.
 *
 * @param[in,out] M     The coefficient matrix.
 * @param[in]     v     The constant vector.
 * @param[out]    x     The variable vector.
 * @param[in]     nEqs  The number of equations.
 *
 * @tparam        T     Type of the elements in `M` and `v`.
 */
template <class T>
void solveLinSystem(T **M, T v[], T x[], const int& nEqs) {
	// Reduce the system to upper triangular form
	gaussElim(M, v, nEqs);

	// Solve the system
	for (int i = nEqs - 1; i >= 0; i--) {
		T temp = v[i];
		for (int j = nEqs - 1; j > i; j--) temp -= x[j] * M[i][j];
		x[i] = temp / M[i][i];
	}
}

/**
 * @brief      Solve a tridiagonal linear system.
 *
 * @param[in]  a     The array with the elements M[i][i - 1] (sub-diagonal).
 * @param[in]  b     The array with the elements M[i][i] (diagonal).
 * @param[in]  c     The array with the elements M[i][i + 1] (sur-diagonal).
 * @param[in]  r     The constant vector.
 * @param[out] x     The variable vector.
 * @param[in]  nEq   The number of equations
 *
 * @tparam     T     Type of the elements in the arrays.
 *
 * @throws     std::invalid_argument  Thrown if nEq > 4096.
 * @throws     std::invalid_argument  Thrown if a[0] isn't NaN.
 * @throws     std::invalid_argument  Thrown if c[-1] isn't NaN.
 */
template <class T>
void tridiagonalSolver(T a[], T b[], T c[], T r[], T x[], const int& nEq) {
	const int maxSize = 4096;
	if (nEq > maxSize) throw std::invalid_argument("Number of equations must not be greater than " + std::to_string(maxSize) + ".");
	if (!isnan(a[0])) throw std::invalid_argument("The first element of a must be nan(\"\").");
	if (!isnan(c[nEq - 1])) throw std::invalid_argument("The last element of c must be nan(\"\").");
	T h[maxSize], p[maxSize];

	h[0] = c[0] / b[0];
	p[0] = r[0] / b[0];

	// Gaussian elimination
	for (int i = 1; i < nEq; i++) {
		T den =  1.0 / (b[i] - a[i] * h[i - 1]);
		h[i] = c[i] * den;
		p[i] = (r[i] - a[i] * p[i - 1]) * den;
	}

	// Backsubsitution
	x[nEq - 1] = p[nEq - 1];
	for (int i = nEq - 2; i >= 0; i--) {
		x[i] = p[i] - h[i] * x[i + 1];
	}
}

/**
 * @brief      Solves a linear boundary value problem in the form y'' = f(x).
 *
 * @param[in]  RHS      The Right Hand Side function (i.e.: f(x)).
 * @param[out] y        The array with the solution of the equation.
 * @param[in]  xL       The integration starting value.
 * @param[in]  xR       The integration stopping value.
 * @param[in]  yL       The left boundary condition.
 * @param[in]  yR       The right boundary condition.
 * @param[in]  nPoints  The number of points to use in the integration.
 *
 * @throws     std::invalid_argument  Thrown if nPoints > 4096 (via
 *                                    tridiagonalSolver()).
 */
void linearBVP(double (*RHS)(const double& x), double y[], const double& xL, const double& xR, const double& yL, const double& yR, const int& nPoints);
