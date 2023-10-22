#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

void std_formula(double, double);

void smart_formula(double, double);

int quadratic_formula1(double, double, double, double&, double&);

int quadratic_formula2(double, double, double, double&, double&);

void print_results(double, double, double, double);

void sortRoots(double&, double&);

void swap(double&, double&);

int sign(double);

int main() {
	cout << setiosflags(ios::scientific);

	double x1 = 1.e12;
	double x2 = 1.e-12;

	sortRoots(x1, x2);

	std_formula(x1, x2);
	cout << endl;
	smart_formula(x1, x2);
	return 0;
}

void std_formula(double x1, double x2) {
	double X1, X2;
	int flag = quadratic_formula1(1, -(x1 + x2), x1 * x2, X1, X2);
	if (flag == 0) {
		cout << "Fools' way" << endl;
		cout << "====================" << endl;

		print_results(x1, x2, X1, X2);
	} else if (flag == 1) {
		cout << "No solutions!" << endl;
	}
}

int quadratic_formula1(double a, double b, double c, double& x1, double& x2) {
	double delta = sqrt(b * b - 4 * a * c);
	if (delta >= 0) {
		x1 = (-b + delta) / (2 * a);
		x2 = (-b - delta) / (2 * a);
		return 0;		
	} else {
		return 1;
	}
}

void smart_formula(double x1, double x2) {
	double X1, X2;
	int flag = quadratic_formula2(1, -(x1 + x2), x1 * x2, X1, X2);
	if (flag == 0) {
		cout << "Smart way" << endl;
		cout << "====================" << endl;
		print_results(x1, x2, X1, X2);
	} else if (flag == 1) {
		cout << "No solutions!" << endl;
	}
}

int quadratic_formula2(double a, double b, double c, double& x1, double& x2) {
	double delta = sqrt(b * b - 4 * a * c);
	if (delta >= 0) {
		x1 = (2 * c) / (-b - sign(b) * delta);
		x2 = (-b - sign(b) * delta) / (2 * a);
		return 0;
	} else {
		return 1;
	}
}

void print_results(double x1, double x2, double X1, double X2) {
	sortRoots(X1, X2);
	cout << "x1 = " << x1 << "\tX1 = " << X1 << "\terr = " << abs(x1 - X1) << endl;
	cout << "x2 = " << x2 << "\tX2 = " << X2 << "\terr = " << abs(x2 - X2) << endl;
}

void sortRoots(double& x1, double& x2) {
	if (x1 > x2) {
		swap(x1, x2);
	}
}

void swap(double& a, double& b) {
	double t = a;
	a = b;
	b = t;
}

int sign(double x) {
	if (x >= 0) return 1;
	else return -1;
}
