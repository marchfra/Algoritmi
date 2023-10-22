#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

// Library functions

double rectangular_quad(double (*)(double), const double, const double, const int);

double midpoint_quad(double (*)(double), const double, const double, const int);

double trapezoidal_quad(double (*)(double), const double, const double, const int);

double simpson_quad(double (*)(double), const double, const double, const int);

double gauss_legendre_quad(double (*)(double), const double, const double, const int = 1, const int = 3);

double horner_pol(const double, const double[], const int);
double horner_pol(const double, const double[], const int, double&);

void swap(double&, double&);
void swap(int&, int&);

// Other functions

double func(double);

int main() {
	cout << setiosflags(ios::scientific);

	ofstream outgnu, outpy;
	outgnu.open("data.dat", ios::out);
	outpy.open("data.csv", ios::out);

	if (!outgnu || !outpy) {
		cout << "Error while opening file." << endl;
		exit(1);
	}

	double x0 = 0.0;
	double dx = 0.1;
	double x = x0;
	double sum = 0.0;

	// x = 0.8;
	// cout << "------------------------------------------------------------" << endl;
	// cout << "h\t\t\t\tM\t\t\t\tT\t\t\t\tG" << endl;
	// cout << "------------------------------------------------------------" << endl;
	// int N = 1;
	// while(N < 10) {
	// 	cout << x/N << "\t" << midpoint_quad(func, 0, x, N) << "\t";
	// 	cout << trapezoidal_quad(func, 0, x, N) << "\t";
	// 	cout << gauss_legendre_quad(func, 0, x, N) << endl;
	// 	N *= 2;
	// }
	// cout << "------------------------------------------------------------" << endl;

	outpy << "x,F" << endl;

	while (x < 25) {
		sum += gauss_legendre_quad(func, x0, x);
		outgnu << x << " " << sum << endl;
		outpy << x << "," << sum << endl;
		x0 = x;
		x += dx;
	}

	outgnu.close();
	outpy.close();

	return 0;
}

double func(double x) {
	if (x != 0.0) return sin(x) / x;
	else return 1;
}

// ==================================================================================================================================================
// Library functions
// ==================================================================================================================================================

double rectangular_quad(double (*F)(double), const double a, const double b, const int n) {
	double sum = 0.0;
	const double h = (b - a) / n;	// Interval width
	double xi = a;

	for (int i = 0; i < n; i++) {	// Rectangular rule
		sum += F(xi);
		xi += h;
	}
	
	return sum * h;					// Finalise the calculation and return
}

double midpoint_quad(double (*F)(double), const double a, const double b, const int n) {
	return gauss_legendre_quad(F, a, b, n, 1);		// Gauss-Legendre rule converges to midpoint rule for Ng = 1
}

double trapezoidal_quad(double (*F)(double), const double a, const double b, const int n) {
	double sum = 0.5 * (F(a) + F(b));	// Take care of the first and last term of the summation
	const double h = (b - a) / n;		// Interval width
	double xi = a;

	for (int i = 1; i < n; i++) {		// Trapezoidal rule
		xi += h;
		sum += F(xi);
	}
	
	return sum * h;						// Finalise the calculation and return
}

double simpson_quad(double (*F)(double), const double a, const double b, const int n) {
	// n must be even!
	if (n % 2 != 0) {
		cout << "simpson_quad(): Invalid argument: n must be even." << endl;
		exit(1);
	}
	
	const double h = (b - a) / n;		// Interval width
	double sum = F(a) + F(b);			// Take care of the first and last term of the summation
	double xi = a;
	double w = 4.0;						// Define the second weight to be 4

	for (int i = 1; i < n; i++) {
		xi += h;
		sum += F(xi) * w;
		w = 6.0 - w;					// Alternate between w = 4 and w = 2
	}
	
	return sum * h / 3;					// Finalise the calculation and return
}

double gauss_legendre_quad(double (*F)(double), const double a, const double b, const int N, const int Ng){
	double weight[8], root[8];							// Later you can try and reduce the memory used
														// Standard C++ requires the static size to be known
														// at compilation time. Fot all we know, Ng could be 
														// a user input, so it can become known only at runtime

	switch(Ng) {
	case 1:
		root[0] = 0;												weight[0] = 2;
		break;
	case 2:
		root[0] =  sqrt(1.0 / 3.0);									weight[0] = 1;
		root[1] = -sqrt(1.0 / 3.0);									weight[1] = 1;
		break;
	case 3:
		root[0] =  0;												weight[0] = 8.0 / 9.0;
		root[1] =  sqrt(3.0 / 5.0);									weight[1] = 5.0 / 9.0;
		root[2] = -sqrt(3.0 / 5.0);									weight[2] = 5.0 / 9.0;
		break;
	case 4:
		root[0] =  sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	weight[0] = (18 + sqrt(30)) / 36;
		root[1] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	weight[1] = (18 + sqrt(30)) / 36;
		root[2] =  sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	weight[2] = (18 - sqrt(30)) / 36;
		root[3] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	weight[3] = (18 - sqrt(30)) / 36;
		break;
	case 5:
		root[0] =  0;												weight[0] = 128.0 / 225.0;
		root[1] =  1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		weight[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0; 
		root[2] = -1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		weight[2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		root[3] =  1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		weight[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		root[4] = -1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		weight[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		break;
	default:
		cout << "gauss_legendre_quad(): Invalid argument: Ng must be in [1, 5]." << endl;
		exit(1);
	}

	const double h = (b - a) / N;						// Interval width
	double sum = 0.0;

	double xc, sumj;

	for (int i = 0; i < N; i++) {						// Loop over intervals
		xc = a + (i + 0.5) * h;							// Define the center of the interval
		sumj = 0.0;										// Initialize sum for this interval
		for (int j = 0; j < Ng; j++) {					// Loop over gaussian points in the interval
			sumj += weight[j] * F(0.5*h*root[j] + xc);	// and apply gaussian rule
		}
		sum += sumj;									// Add partial sum to total integral
	}

	return 0.5 * h * sum;								// Finalize calculation and return
}

double horner_pol(const double x, const double a[], const int degree) {
	double p = a[degree];
	for (int i = degree-1; i >= 0; i--) {
		p = p * x + a[i];
	}
	return p;
}

double horner_pol(const double x, const double a[], const int degree, double& dpdx) {
	double p = a[degree];
	dpdx = 0.0;
	for (int i = degree-1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p = p * x + a[i];
	}
	return p;
}

void swap(double& a, double& b) {
	double t = a;
	a = b;
	b = t;
}

void swap(int& a, int& b) {
	int t = a;
	a = b;
	b = t;
}
