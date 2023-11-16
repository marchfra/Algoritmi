/**
 * @file swap.hpp
 * 
 * @brief Implementation of the swap utility function.
 */

#pragma once

/**
 * @brief Swap two values.
 * 
 * This function swaps the two values in input.
 * 
 * @tparam a, b The values to be swapped.
 */
template <class T>
void swap(T& a, T& b) {
	T temp = b;
	b = a;
	a = temp;
}
