#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

// Library functions

double rectangular_quad(double (*)(double), const double, const double, const int);

double midpoint_quad(double (*)(double), const double, const double, const int);

double trapezoidal_quad(double (*)(double), const double, const double, const int);

double simpson_quad(double (*)(double), const double, const double, const int);

double gauss_legendre_quad(double (*)(double), const double, const double, const int = 1, const int = 3);
// Integrate over a rectangle
double gauss_legendre_quad2D(double (*)(double, double), const double, const double, const double, const double, const int = 1, const int = 3);

double horner_pol(const double, const double[], const int);
double horner_pol(const double, const double[], const int, double&);

void swap(double&, double&);
void swap(int&, int&);

void time_test(void (*)(), const string, const long int = 1e6);
void time_test(double (*)(), const string, const long int = 1e6);
void time_test_monte_carlo(double (*)(int), const int, const string, const long int = 1e6);
// Time test for Rectangular, Midpoint, Trapezoid, Simpson quad
void time_test(double (*)(double (*)(double), double, double, int), double(*)(double), const double, const double, int, const string, const long int = 1e6);
// Time test for Gauss quad
void time_test(double (*)(double (*)(double), double, double, int, int), double(*)(double), const double, const double, int, int, const string, const long int = 1e6);


// Other functions

double pi_monte_carlo(int);

int main() {
	srand48(time(NULL));

	ofstream out;
	out.open("data.csv");

	cout.precision(8);
	cout << "π_esp = " << pi_monte_carlo(1e6) << endl;
	cout << "π =     " << M_PI << endl << endl;

	int N = 2;

	double err = 1.0;
	double tol = 1.0e-4;
	double I;

	while (err > tol) {
		I = pi_monte_carlo(N);
		err = fabs(I / M_PI - 1.0);
		N++;
	}
	cout << "N such that err < " << tol << ": " << N - 1 << endl;

	out << "N,err,sigma,pi" << endl;

	N = 4;

	while (N < (int)1e7) {
		I = pi_monte_carlo(N);
		err = fabs(I / M_PI - 1.0);
		double var = (double)N / (N - 1) * (4.0 * I - I * I);
		double sigma = sqrt(var / N);
		out << N << "," << err << "," << sigma << "," << I << endl;
		N *= 2;
	}


	out.close();

	return 0;
}

double pi_monte_carlo(int N) {
	double x, y;
	int counter = 0;
	for (int i = 0; i < N; i++) {
		x = -1.0 + 2.0 * drand48();  y = -1.0 + 2.0 * drand48();
		if (x*x + y*y <= 1) counter++;
	}

	return 4.0 * (double)counter / N;
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

double gauss_legendre_quad(double (*F)(double), const double a, const double b, const int N, const int Ng) {
	double weight[8], root[8];

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
		xc = a + (i + 0.5) * h;							// Define the centre of the interval
		sumj = 0.0;										// Initialize sum for this interval
		for (int j = 0; j < Ng; j++) {					// Loop over gaussian points in the interval
			sumj += weight[j] * F(0.5*h*root[j] + xc);	// and apply gaussian rule
		}
		sum += sumj;									// Add partial sum to total integral
	}

	return 0.5 * h * sum;								// Finalize calculation and return
}

double gauss_legendre_quad2D(double (*F)(double, double), const double xa, const double xb, const double ya, const double yb, const int N, const int Ng) {
	double weight[8], root[8];

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

	const double xh = (xb - xa) / N;					// Interval width along x
	const double yh = (yb - ya) / N;					// Interval width along y
	double sum = 0.0;

	double xc, yc, xsum, ysum;

	for (int i = 0; i < N; i++) {						// Loop over x intervals
		for (int j = 0; j < N; j++) {					// Loop over y intervals
			xc = xa + (i + 0.5) * xh;					// Define the centre of the
			yc = ya + (j + 0.5) * yh;					// x and y intervals

			ysum = 0.0;									// Initialize sum over x for this interval
			for (int jk = 0; jk < Ng; jk++) {			// Loop over y gaussian points
				xsum = 0.0;								// Initialize sum over y for this interval
				for (int ik = 0; ik < Ng; ik++) {		// Loop over x gaussian points and apply gaussian rule
					xsum += weight[ik] * F(0.5*xh*root[ik] + xc, 0.5*yh*root[jk] + yc);
				}
				ysum += weight[jk] * xsum;				// Add weighted xsum to y sum
			}
			sum += ysum;								// Add partial sum to total integral
		}
	}

	return 0.25 * xh * yh * sum;						// Finalise calculation and return
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

void time_test(void (*f)(), const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(), const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(int), const int N, const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(N);
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(double (*)(double), double, double, int), double(*func)(double), const double a, const double b, int N, const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N);
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(double (*)(double), double, double, int, int), double(*func)(double), const double a, const double b, int N, int Ng, const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N, Ng);
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}
