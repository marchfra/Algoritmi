#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

double f1(double x) {
	return x*x*x + 5.0*x*x - 3.0*x + 1.0;
}

double f2(double x) {
	return 2.0*x*x*x - 6.0*x*x + 2.0*x - 1.0;
}

double f3(double x) {
	return 2.0*x*x*x + 3.0*x + 1.0;
}

double horner_pol(const double, const double[], const int);

int main() {
	double x = -4.0;
	double dx = 0.1;

	#if 0
		double a[] = {1.0, 5.0, -3.0, 1.0};
		int n = sizeof(a) / sizeof(a[0]) - 1;
		cout << "x,pol,horn" << endl;
		while (x <= 4.0) {
			cout << x << "," << f1(x) << "," << horner_pol(x, a, n) << endl;
			x += dx;
		}
	#endif

	#if 0
		double a[] = {-1.0, 2.0, -6.0, 2.0};
		int n = sizeof(a) / sizeof(a[0]) - 1;
		cout << "x,pol,horn" << endl;
		while (x <= 4.0) {
			cout << x << "," << f2(x) << "," << horner_pol(x, a, n) << endl;
			x += dx;
		}
	#endif

	#if 1
		double a[] = {1.0, 3.0, 0, 2.0};
		int n = sizeof(a) / sizeof(a[0]) - 1;
		cout << "x,pol,horn" << endl;
		while (x <= 4.0) {
			cout << x << "," << f2(x) << "," << horner_pol(x, a, n) << endl;
			x += dx;
		}
	#endif

	return 0;
}

double horner_pol(const double x, const double a[], const int degree) {
	double p = a[degree];
	for (int i = degree-1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}
