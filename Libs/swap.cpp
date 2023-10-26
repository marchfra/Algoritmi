#include "debug.hpp"
#include "swap.hpp"

void swap(int& a, int& b) {
	int temp = a;
	a = b;
	b = temp;
}

void swap(float& a, float& b) {
	float temp = a;
	a = b;
	b = temp;
}

void swap(double& a, double& b) {
	double temp = a;
	a = b;
	b = temp;
}
