// Francesco Marchisotti (fila C)
// Date: 02/11/2023

// Code output:
// ******************************************************************
// Gaussian:    n =   256; Wx = 471.827; Wy = -20.359; Wtot = 451.468
// Trapezoidal: n = 32768; Wx = 471.827; Wy = -20.359; Wtot = 451.468
// ******************************************************************


#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"

using std::cout;
using std::endl;

double IntegralX(double dTheta);
double IntegralY(double dTheta);

double Integrand(double dTheta);


int main() {
	double dThetaMin = 0.0;
	double dThetaMax = 4.0 * M_PI;
	int nN = 4;
	int nNg = 3;

	double err = std::numeric_limits<double>::max();
	double tol = 1.0e-6;
	double dWGauss, dWGauss_ = 0.0;
	double dWXGauss, dWYGauss;
	// cout << "---------- Guassian Integral ----------" << endl;
	while (err > tol) {
		dWXGauss = gauss_legendre_quad(IntegralX, dThetaMin, dThetaMax, nN, nNg);
		dWYGauss = gauss_legendre_quad(IntegralY, dThetaMin, dThetaMax, nN, nNg);
		dWGauss = dWXGauss + dWYGauss;
		err = fabs(dWGauss - dWGauss_);	
		dWGauss_ = dWGauss;
		nN *= 2;
		// cout << "N = " << std::setw(5) << nN << ";   W = " << dWGauss << ";   err = " << err << endl;	
	}

	cout << "Gaussian:    n =   " << nN << "; Wx = " << dWXGauss << "; Wy = " << dWYGauss << "; Wtot = " << dWGauss << endl;

	nN = 4;
	err = std::numeric_limits<double>::max();
	double dWTrap, dWTrap_ = 0.0;
	double dWXTrap, dWYTrap;
	// cout << "\n-------- Trapezoidal Integral --------" << endl;
	while (err > tol) {
		dWXTrap = trapezoidal_quad(IntegralX, dThetaMin, dThetaMax, nN);
		dWYTrap = trapezoidal_quad(IntegralY, dThetaMin, dThetaMax, nN);
		dWTrap = dWXTrap + dWYTrap;
		err = fabs(dWTrap - dWTrap_);	
		dWTrap_ = dWTrap;
		nN *= 2;
		// cout << "N = " << std::setw(5) << nN << ";   W = " << dWTrap << ";   err = " << err << endl;	
	}

	cout << "Trapezoidal: n = " << nN << "; Wx = " << dWXTrap << "; Wy = " << dWYTrap << "; Wtot = " << dWTrap << endl;

	cout << endl << "Actual integral" << endl;
	nN = 4;
	err = std::numeric_limits<double>::max();
	dWTrap_ = 0.0;
	// cout << "\n-------- Trapezoidal Integral --------" << endl;
	while (err > tol) {
		dWTrap = trapezoidal_quad(Integrand, dThetaMin, dThetaMax, nN);
		err = fabs(dWTrap - dWTrap_);
		dWTrap_ = dWTrap;
		nN *= 2;
		// cout << "N = " << std::setw(5) << nN << ";   W = " << dWTrap << ";   err = " << err << endl;	
	}

	cout << "Trapezoidal: n = " << nN << "; Wtot = " << dWTrap << endl;

	nN = 4;
	err = std::numeric_limits<double>::max();
	dWGauss_ = 0.0;
	// cout << "---------- Guassian Integral ----------" << endl;
	while (err > tol) {
		dWGauss = gauss_legendre_quad(Integrand, dThetaMin, dThetaMax, nN, nNg);
		err = fabs(dWGauss - dWGauss_);	
		dWGauss_ = dWGauss;
		nN *= 2;
		// cout << "N = " << std::setw(5) << nN << ";   W = " << dWGauss << ";   err = " << err << endl;	
	}

	cout << "Gaussian:    n =   " << nN << "; Wtot = " << dWGauss << endl;
	cout << "Exact = " << sqrt(17) - 1 << endl;

	return 0;
}

double IntegralX(double dTheta) {
	double dB = M_PI;
	double dArg = dTheta / dB;
	double dRoot = sqrt(1 + dArg*dArg*dArg*dArg);

	double dTerm1 = 2.0 / dB*dB*dB*dB * dTheta*dTheta*dTheta * 1.0 / (dRoot*dRoot*dRoot) * cos(dTheta)*cos(dTheta);
	double dTerm2 = 1.0 / dRoot * sin(dTheta) * cos(dTheta);

	return dTerm1 + dTerm2;
}

double IntegralY(double dTheta) {
	double dB = M_PI;
	double dArg = dTheta / dB;
	double dRoot = sqrt(1 + dArg*dArg*dArg*dArg);

	double dTerm1 = 2.0 / dB*dB*dB*dB * dTheta*dTheta*dTheta * 1.0 / (dRoot*dRoot*dRoot) * sin(dTheta) * cos(dTheta);
	double dTerm2 = -1.0 / dRoot * cos(dTheta)*cos(dTheta);

	return dTerm1 + dTerm2;
}

double Integrand(double dTheta) {
	double dB4 = M_PI*M_PI*M_PI*M_PI;
	double dTheta3 = dTheta*dTheta*dTheta;
	return 2.0 / dB4 * dTheta3 / sqrt(1 + dTheta*dTheta3);
}

// =====================================================================================================================
// Library Functions
// =====================================================================================================================
/*
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
		std::cout << "gauss_legendre_quad(): Invalid argument: Ng must be in [1, 5]." << std::endl;
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
*/
