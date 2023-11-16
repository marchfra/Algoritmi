#include "swap.hpp"

template <class T>
int partition(T arr[], const int& first, const int& last) {
	srand(time(NULL));
	// int pivot = rand() % last;
	T x = arr[last];
	int i = first - 1;

	for (int j = first; j < last; j++) {
		if (arr[j] <= x) {
			i++;
			swap(arr[i], arr[j]);
		}
	}
	swap(arr[i + 1], arr[last]);
	return i + 1;
}

template <class T>
void quicksort(T arr[], const int& size, const int m = -1) {

	// //////////////////////////////////////////////
	//
	// Sort an array using quicksort method
	// 
	// arr   [in-out] : the array to sort
	// size  [in]     : the size of the array
	// m     [in]     : the number of elements to sort
	//                 if omitted sort whole array
	// 
	// Last modified: 07 Nov 2023
	//
	// //////////////////////////////////////////////
}
