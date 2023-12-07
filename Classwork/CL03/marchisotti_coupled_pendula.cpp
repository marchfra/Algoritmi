/*
 * Name: Francesco Marchisotti
 * Date: 07/12/2023
 *
 * Code output:
 * Loop break at nStep =   1047; t =     104.7
 *               err1 =  0.2504; err2 = 0.0293
 *               tp1 =      101; tp2 =     103
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"

#define MAX(a, b) (a > b ? a : b)

using std::cout;
using std::endl;

const double g_L = 1.0;
const double g_G = 9.8;
const double g_K = 0.8;

void RHS(const double& t, double Y[], double R[]);
void RHS(double Y[], double R[]);

int main() {
	const double t0 = 0.0;
	const double dt = 0.1;

	// Initial conditions
	const double theta10 = -0.25 * M_PI;
	const double theta20 = 0.0;
	const double omega10 = 0.0;
	const double omega20 = 0.0;
	const double y0[] = {theta10, theta20, omega10, omega20};

	const int nEq = static_cast<int>(sizeof(y0) / sizeof(y0[0]));
	double yRK[nEq];
	double yPV[nEq];
	for (int i = 0; i < nEq; i++) {
		yRK[i] = y0[i];
		yPV[i] = y0[i];
	}

	const double tol = 0.25;
	double err1 = std::numeric_limits<double>::min();
	double err2 = std::numeric_limits<double>::min();

	// Number of times the pendulum changes direction
	int nTurningPoints1 = 0;
	int nTurningPoints2 = 0;
	double prevTheta1 = theta10;
	double prevTheta2 = theta20;

	std::ofstream out;
	out.open("../data/data.csv");
	if (!out) exit(5);

	out << "t,th1RK,th2RK,w1RK,w2RK,th1PV,th2PV,w1PV,w2PV" << endl;

	double t = t0;
	int nStep = 0;
	out << t << "," << yRK[0] << "," << yRK[1] << "," << yRK[2] << "," << yRK[3]
             << "," << yPV[0] << "," << yPV[1] << "," << yPV[2] << "," << yPV[3] << endl;

	while (MAX(err1, err2) <= tol) {
		// Advance solutions
		rk4Step(t, yRK, RHS, dt, nEq);
		pVerlet(t, yPV, RHS, dt, nEq);
		t += dt;

		out << t << "," << yRK[0] << "," << yRK[1] << "," << yRK[2] << "," << yRK[3]
		         << "," << yPV[0] << "," << yPV[1] << "," << yPV[2] << "," << yPV[3] << endl;

		// Check if either pendulum changed direction
		if (prevTheta1 * yRK[0] <= 0.0)
			nTurningPoints1++;
		if (prevTheta2 * yRK[1] <= 0.0)
			nTurningPoints2++;

		// Update previous angles
		prevTheta1 = yRK[0];
		prevTheta2 = yRK[1];

		// Compute error
		err1 = fabs(yRK[0] - yPV[0]) / M_PI;
		err2 = fabs(yRK[1] - yPV[1]) / M_PI;

		nStep++;
	}

	out.close();

	cout.precision(4);
	cout << "Loop break at nStep = " << std::setw(6) << nStep           << "; t =    " << std::setw(6) << t << endl
	     << "              err1 =  " << std::setw(6) << err1            << "; err2 = " << std::setw(6) << std::setprecision(3) << err2 << endl
	     << "              tp1 =   " << std::setw(6) << nTurningPoints1 << "; tp2 =  " << std::setw(6) << nTurningPoints2 << endl;
	return 0;
}

void RHS(const double& t, double Y[], double R[]) {
	RHS(Y, R);
}

void RHS(double Y[], double R[]) {
	double theta1 = Y[0];
	double theta2 = Y[1];
	double omega1 = Y[2];
	double omega2 = Y[3];

	// dtheta = omega
	// domega = f
	R[0] = omega1;
	R[1] = omega2;
	R[2] = -g_G / g_L * sin(theta1) - g_K * (sin(theta1) - sin(theta2));
	R[3] = -g_G / g_L * sin(theta2) + g_K * (sin(theta1) - sin(theta2));
}

// void rk4Step(const double &t, double Y[],
//              void (*RHSFunc)(const double &t, double Y[], double RHS[]),
//              const double &dt, const int &neq) {
//   if (neq > 64)
//     throw std::invalid_argument("neq must be less than 64");

//   double Ystar[64], k1[64], k2[64], k3[64], k4[64];

//   RHSFunc(t, Y, k1);

//   for (int i = 0; i < neq; i++) {
//     Ystar[i] = Y[i] + 0.5 * dt * k1[i];
//   }
//   RHSFunc(t + 0.5 * dt, Ystar, k2);

//   for (int i = 0; i < neq; i++) {
//     Ystar[i] = Y[i] + 0.5 * dt * k2[i];
//   }
//   RHSFunc(t + 0.5 * dt, Ystar, k3);

//   for (int i = 0; i < neq; i++) {
//     Ystar[i] = Y[i] + dt * k3[i];
//   }
//   RHSFunc(t + dt, Ystar, k4);

//   for (int i = 0; i < neq; i++) {
//     Y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
//   }
// }

// void pVerlet(const double &t, double Y[],
//              void (*RHSFunc)(double Y[], double RHS[]), const double &dt,
//              const int &nEq) {
//   const int nMax = 64;
//   if (nEq > 64)
//     throw std::invalid_argument("nEq must be less than " +
//                                 std::to_string(nMax) + ".");
//   if (nEq % 2 != 0)
//     throw std::invalid_argument("nEq must be even");

//   double a[nMax / 2];

//   int nParticles = nEq / 2;

//   double *x = Y;
//   double *v = Y + nParticles;

//   for (int i = 0; i < nParticles; i++) {
//     x[i] += 0.5 * dt * v[i];
//   }

//   RHSFunc(Y, a);

//   for (int i = 0; i < nParticles; i++) {
//     v[i] +=
//         dt *
//         (a + nParticles)[i]; // Note: index of a is the index of the velocities
//     x[i] += 0.5 * dt * v[i];
//   }
// }
