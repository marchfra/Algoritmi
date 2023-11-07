#include <iostream>

#include "../include/debug.hpp"
#include "../include/derivative.hpp"

void EulerStep(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");

	double rhs[64];

	RHSFunc(t, Y, rhs);
	for (int i = 0; i < neq; i++) {
		Y[i] += dt * rhs[i];
	}
}

void RK2Step(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq, bool rule) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");
	double Ystar[64], k1[64], k2[64];

	if (!rule) {
		// Midpoint
		RHSFunc(t, Y, k1);

		for (int i = 0; i < neq; i++) {
			Ystar[i] = Y[i] + 0.5 * dt * k1[i];
		}
		RHSFunc(t + 0.5 * dt, Ystar, k2);

		for (int i = 0; i < neq; i++) {
			Y[i] += dt * k2[i];
		}
	} else {
		// Modified Euler
		RHSFunc(t, Y, k1);

		for (int i = 0; i < neq; i++) {
			Ystar[i] = Y[i] + dt * k1[i];
		}
		RHSFunc(t + dt, Ystar, k2);

		for (int i = 0; i < neq; i++) {
			Y[i] += dt * 0.5 * (k1[i] + k2[i]);
		}
	}
}

void RK4Step(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");

	double Ystar[64], k1[64], k2[64], k3[64], k4[64];

	RHSFunc(t, Y, k1);

	for (int i = 0; i < neq; i++) {
		Ystar[i] = Y[i] + 0.5 * dt * k1[i];
	}
	RHSFunc(t + 0.5 * dt, Ystar, k2);

	for (int i = 0; i < neq; i++) {
		Ystar[i] = Y[i] + 0.5 * dt * k2[i];
	}
	RHSFunc(t + 0.5 * dt, Ystar, k3);

	for (int i = 0; i < neq; i++) {
		Ystar[i] = Y[i] + dt * k3[i];
	}
	RHSFunc(t + dt, Ystar, k4);

	for (int i = 0; i < neq; i++) {
		Y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
	}
}

