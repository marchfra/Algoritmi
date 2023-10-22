#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

void example1();
float f11(float);
float f12(float);
float f13(float);

void example2();
float f21(float);
float f22(float);
float f23(float);


int main() {
	cout << setiosflags(ios::scientific);

	example1();

	cout << endl;

	example2();

	return 0;
}

void example1() {
	cout << "Example #1: compute sqrt(x^2 + 1) - x for large x" << endl;
	cout << "=================================================" << endl;
	int size = 7;
	float x = 1.e4;

	for (int i = 0; i < size; i++) {
		cout << "x = " << x << ";\tfx1 = " << f11(x) << ";\tfx2 = " << f12(x) << ";\tfx3 = " << f13(x) << endl;
		x *= 10;
	}
}

float f11(float x) {
	return sqrt(x * x + (float)1) - x;
}

float f12(float x) {
	float num = 1;
	float den = sqrt(x * x + (float)1) + x;
	return num / den;
}

float f13(float x) {
	return 1 / (2*x) - 1 / (8*x*x*x);
}

// ============================================================

void example2() {
	cout << "Example #2: compute 1-cos(x) for small x" << endl;
	cout << "========================================" << endl;
	int size = 9;
	float x = 1.e-1;

	for (int i = 0; i < size; i++) {
		cout << "x = " << x << ";\tfx1 = " << f21(x) << ";\tfx2 = " << f22(x) << ";\tfx3 = " << f23(x) << endl;
		x *= 0.1;
	}
}

float f21(float x) {
	return 1 - cos(x);
}

float f22(float x) {
	float num = sin(x) * sin(x);
	float den = 1 + cos(x);
	return num / den;
}

float f23(float x) {
	return x*x / 2 - x*x*x*x / 24;
}
