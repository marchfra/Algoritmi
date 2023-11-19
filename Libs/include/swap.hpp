/**
 * @file swap.hpp
 *
 * @brief      Implementation of the swap utility function.
 *
 * @author     Francesco Marchisotti
 *
 * @date       19/11/2023
 */
#pragma once

/**
 * @brief      Swap two values.
 *
 * This function swaps the two values in input.
 *
 * @param      a     The first value.
 * @param      b     The second value.
 *
 * @tparam     T     Type of elements to be swapped.
 */
template <class T>
void swap(T& a, T& b) {
	T temp = b;
	b = a;
	a = temp;
}
