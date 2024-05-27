/**
 * @file    lin_alg.hpp
 *
 * @brief   This file implements linear algorithms.
 *
 * @author  Francesco Marchisotti
 *
 * @date    2023-11-24
 */

#pragma once

#include "../include/swap.hpp"

#include <iomanip>
#include <iostream>

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
