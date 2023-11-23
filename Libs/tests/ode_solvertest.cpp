#include <cmath>

#include "test_config.hpp"
#include "../include/ode_solver.hpp"

void RHS1(const double& t, double Y[], double R[]);

double exact1(const double& t);

TEST_CASE("testing eulerStep function") {
	
}

double exact1(const double& t) {
	return exp(-0.5 * t*t);
}

void RHS1(const double& t, double Y[], double R[]) {
	// dy/dt = -ty
	R[0] = -t * Y[0];
}
