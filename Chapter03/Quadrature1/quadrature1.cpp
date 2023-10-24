#include <iostream>
#include <iomanip>
#include <cmath>

// Commenta bene tutto

using namespace std;

double rectangular_quad(double (*)(double), double, double, int);

double trapezoidal_quad(double (*)(double), double, double, int);

double simpson_quad(double (*)(double), double, double, int);

double gauss_legendre_quad(double (*)(double), double, double, int, int);

void convergence_test(double (*)(double (*)(double), double, double, int), double(*)(double), double, double, double, int);

double func(double);

double exact(double, double);

int main() {

	cout << setiosflags(ios::scientific);

	double a = 0, b = 1;
	int N = 4;

	double rect = rectangular_quad(func, a, b, N);
	double trap = trapezoidal_quad(func, a, b, N);
	double simp = simpson_quad(func, a, b, N);

	cout << "Exact:       " << exact(a, b) << endl;
	cout << "Rectangular: " << rect << "; (N = " << N << ")" << endl;
	cout << "Trapezoidal: " << trap << "; (N = " << N << ")" << endl;
	cout << "Simpson:     " << simp << "; (N = " << N << ")" << endl;

	cout << endl;

	cout << "Convergence test" << endl;
	cout << "================" << endl;
	double tol = 1.0e-5;
	cout << "Rectangular: ";
	convergence_test(rectangular_quad, func, a, b, tol, N);
	cout << "Trapezoidal: ";
	convergence_test(trapezoidal_quad, func, a, b, tol, N);
	cout << "Simpson:     ";
	convergence_test(simpson_quad, func, a, b, tol, N);

	return 0;
}

void convergence_test(double (*rule)(double (*)(double), double, double, int), double (*F)(double), double a, double b, double tol, int N) {
	double err = 1.0;
	double I_N = rule(F, a, b, N), I_NN;
	while (err >= tol) {
		N *= 2;
		I_NN = rule(F, a, b, N);
		err = fabs(I_NN - I_N);
		I_N = I_NN;

		if (N >= 2000000) {
			cout << "Error: N > 2e6; limit exceeded)" << endl;
		}
	}

	cout << I_NN << "; N = " << N << endl;
}

double rectangular_quad(double (*F)(double), double a, double b, int n) {
	double sum = 0.0;
	double h = (b - a) / n;
	double xi = a;

	for (int i = 0; i < n; i++) {
		sum += F(xi);
		xi += h;
	}
	
	return sum * h;
}

double trapezoidal_quad(double (*F)(double), double a, double b, int n) {
	double sum = 0.5 * (F(a) + F(b));	// Take care of the first and last term of the summation
	double h = (b - a) / n;
	double xi = a;

	for (int i = 1; i < n; i++) {
		xi += h;
		sum += F(xi);
	}
	
	return sum * h;						// Finalise the calculation and return
}

double simpson_quad(double (*F)(double), double a, double b, int n) {
	// n must be even!
	if (n % 2 != 0) {
		throw invalid_argument("n must be even!");
	}
	
	double sum = F(a) + F(b);			// Take care of the first and last term of the summation
	double h = (b - a) / n;
	double xi = a;
	double w = 4.0;						// Define the second weight to be 4

	for (int i = 1; i < n; i++) {
		xi += h;
		sum += F(xi) * w;
		w = 6.0 - w;					// Alternate between w = 4 and w = 2
	}
	
	return sum * h / 3;					// Finalise the calculation and return
}

double func(double x) {
	return exp(-x);
}

double exact(double a, double b) {
	return func(a) - func(b);
}
