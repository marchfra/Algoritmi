#include <iostream>

#include "../include/mycode.hpp"

int factorial(int n) {
	if (n < 0) {
		throw std::invalid_argument("n must be non-negative");
	} else if (n == 0) {
		return 1;
	} else {
		return n * factorial(n - 1);
	}
}

bool returnTrue(int flag) {
	return flag == 0 ? true : false;
}
