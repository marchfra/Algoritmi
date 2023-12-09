#include "../include/ode_solver.hpp"

#include "../include/debug.hpp"

void eulerStep(const double &t, double Y[],
               void (*RHSFunc)(const double &t, double Y[], double RHS[]),
               const double &dt, const int &neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");

	double rhs[64];

	RHSFunc(t, Y, rhs);
	for (int i = 0; i < neq; i++) Y[i] += dt * rhs[i];
}

void rk2Step(const double &t, double Y[],
             void (*RHSFunc)(const double &t, double Y[], double RHS[]),
             const double &dt, const int &neq, bool rule) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");
	double Ystar[64], k1[64], k2[64];

	if (!rule) {
		// Midpoint
		RHSFunc(t, Y, k1);

		for (int i = 0; i < neq; i++) Ystar[i] = Y[i] + 0.5 * dt * k1[i];
		RHSFunc(t + 0.5 * dt, Ystar, k2);

		for (int i = 0; i < neq; i++) Y[i] += dt * k2[i];
	} else {
		// Modified Euler
		RHSFunc(t, Y, k1);

		for (int i = 0; i < neq; i++) Ystar[i] = Y[i] + dt * k1[i];
		RHSFunc(t + dt, Ystar, k2);

		for (int i = 0; i < neq; i++) Y[i] += dt * 0.5 * (k1[i] + k2[i]);
	}
}

void rk4Step(const double &t, double Y[],
             void (*RHSFunc)(const double &t, double Y[], double RHS[]),
             const double &dt, const int &neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");

	double Ystar[64], k1[64], k2[64], k3[64], k4[64];

	RHSFunc(t, Y, k1);

	for (int i = 0; i < neq; i++) Ystar[i] = Y[i] + 0.5 * dt * k1[i];
	RHSFunc(t + 0.5 * dt, Ystar, k2);

	for (int i = 0; i < neq; i++) Ystar[i] = Y[i] + 0.5 * dt * k2[i];
	RHSFunc(t + 0.5 * dt, Ystar, k3);

	for (int i = 0; i < neq; i++) Ystar[i] = Y[i] + dt * k3[i];
	RHSFunc(t + dt, Ystar, k4);

	for (int i = 0; i < neq; i++)
		Y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

void pVerlet(const double &t, double Y[],
             void (*RHSFunc)(double Y[], double RHS[]), const double &dt,
             const int &nEq) {
	const int nMax = 64;
	if (nEq > 64)
		throw std::invalid_argument("nEq must be less than " +
		                            std::to_string(nMax) + ".");
	if (nEq % 2 != 0) throw std::invalid_argument("nEq must be even");

	double a[nMax / 2];

	int nParticles = nEq / 2;

	double *x = Y;
	double *v = Y + nParticles;

	for (int i = 0; i < nParticles; i++) x[i] += 0.5 * dt * v[i];

	RHSFunc(Y, a);

	for (int i = 0; i < nParticles; i++) {
		v[i] +=
			dt *
			(a +
		     nParticles)[i];  // Note: index of a is the index of the velocities
		x[i] += 0.5 * dt * v[i];
	}
}

void vVerlet(const double &t, double Y[],
             void (*RHSFunc)(double Y[], double RHS[]), const double &dt,
             const int &nEq) {
	const int nMax = 64;
	if (nEq > 64)
		throw std::invalid_argument("nEq must be less than " +
		                            std::to_string(nMax) + ".");
	if (nEq % 2 != 0) throw std::invalid_argument("nEq must be even");

	double a[nMax / 2];

	int nParticles = nEq / 2;

	double *x = Y;
	double *v = Y + nParticles;

	RHSFunc(Y, a);

	for (int i = 0; i < nParticles; i++) {
		v[i] +=
			0.5 * dt *
			(a +
		     nParticles)[i];  // Note: index of a is the index of the velocities
		x[i] += dt * v[i];
	}

	RHSFunc(Y, a);

	for (int i = 0; i < nParticles; i++) {
		v[i] +=
			0.5 * dt *
			(a +
		     nParticles)[i];  // Note: index of a is the index of the velocities
	}
}
