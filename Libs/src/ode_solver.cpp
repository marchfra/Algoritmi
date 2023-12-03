#include "../include/ode_solver.hpp"

#include "../include/debug.hpp"

void eulerStep(const double& t, double Y[],
               void (*RHSFunc)(const double& t, double Y[], double RHS[]),
               const double& dt, const int& neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");

	double rhs[64];

	RHSFunc(t, Y, rhs);
	for (int i = 0; i < neq; i++) {
		Y[i] += dt * rhs[i];
	}
}

void rk2Step(const double& t, double Y[],
             void (*RHSFunc)(const double& t, double Y[], double RHS[]),
             const double& dt, const int& neq, bool rule) {
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

void rk4Step(const double& t, double Y[],
             void (*RHSFunc)(const double& t, double Y[], double RHS[]),
             const double& dt, const int& neq) {
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

void pVerlet(const double& t, double Y[],
             void (*RHSFunc)(double Y[], double RHS[]), const double& dt,
             const int& neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");
	if (neq % 2 != 0) throw std::invalid_argument("neq must be even");

	double a[32];

	// int npos = neq / 2;

	// double *x = Y;
	// double *v = Y + npos;

	// for (int i = 0; i < npos; i++) {
	// 	x[i] += 0.5 * dt * v[i];
	// }

	// RHSFunc(Y, a);

	// for (int i = 0; i < npos; i++) {
	// 	v[i] += dt * a[i];
	// 	x[i] += 0.5 * dt * v[i];
	// }


	int x_start = 0, x_stop = neq / 2;
	int v_start = x_stop, v_stop = neq;

	for (int i = x_start; i < x_stop; i++) {
		Y[i] += 0.5 * dt * Y[i + v_start];
	}

	RHSFunc(Y, a);

	for (int i = v_start; i < v_stop; i++) {
		Y[i] += dt * a[i];
		Y[i - v_start] += 0.5 * dt * Y[i];
	}
}

void vVerlet(const double& t, double Y[],
             void (*RHSFunc)(double Y[], double RHS[]), const double& dt,
             const int& neq) {
	if (neq > 64) throw std::invalid_argument("neq must be less than 64");
	if (neq % 2 != 0) throw std::invalid_argument("neq must be even");

	int x_start = 0, x_stop = neq / 2;
	int v_start = x_stop, v_stop = neq;

	double a[64];

	RHSFunc(Y, a);

	for (int i = v_start; i < v_stop; i++) {
		Y[i] += 0.5 * dt * a[i - v_start];
		Y[i - v_start] += dt * Y[i];
	}

	RHSFunc(Y, a);

	for (int i = v_start; i < v_stop; i++) {
		Y[i] += 0.5 * dt * a[i - v_start];
	}
}
