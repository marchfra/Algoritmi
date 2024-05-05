/**
 * @file main.cpp
 *
 * @author Francesco Marchisotti
 *
 * @brief Main file for the course's final project
 *
 * @date 2024-05-04
 */

#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/lin_alg.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/ode_solver.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

// Problem data
const static double B  = 4.0e-5;  //!< Drag coefficient
const static double V0 = 10.0;    //!< Initial velocity
const static double L  = 10.0;    //!< Target distance
double gTheta;                    //!< Initial launch angle

// Dimensional factors
const static double chi = L;
const static double mu  = 1.0;
const static double g   = 9.81;
const static double tau = sqrt(chi / g);
const static double b   = B * chi / mu;
const static double v0  = V0 * tau / chi;

/**
 * @brief Right Hand Side of the ODEs system
 *
 * @param[in] Y The input values of the variables
 * @param[out] R The output values of the variables
 */
void RHS(double Y[], double R[]);

/**
 * @brief Residual function of x(t)
 *
 * This function takes in input a time value and integrates the system of ODEs
 * upto (and including) that time, and returns the value of x(t) - 1 in order to
 * find the time at which x == 1.
 *
 * @param[in] time The time at which to evaluate x
 *
 * @return x(t) - 1
 */
double x_t(const double& time);

/**
 * @brief Residual function of y(t)
 *
 * This function takes in input a time value and integrates the system of ODEs
 * upto (and including) that time, and returns the value of y(t*) in order to
 * find the x-position at which y == 0.
 *
 * @param[in] time The time at which to evaluate y
 *
 * @return y(t)
 */
double y_t(const double& time);

int main() {
	std::ofstream out;

	// Constants output
	out.open("data/constants.csv");
	if (!out) exit(5);

	out << "chi,tau,mu,B,b" << endl;
	out << chi << "," << tau << "," << mu << "," << B << "," << b << endl;

	out.close();


	// Integration range
	const double t0    = 0.0;   // Starting point
	const double t_max = 10.0;  // Ending point
	const int nStep    = 1000;  // Number of steps
	const double dt    = (t_max - t0) / nStep;

	// Initialisation
	double t          = t0;
	const double y0[] = {0.0, 0.0, v0 * cos(gTheta), v0 * sin(gTheta)};
	const int nEq =
		static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));

	// Initial angle
	const double theta0    = 0.68;  // Initial guess for launch angle
	const double theta_max = 0.9;   // Max value for theta
	const int nTheta       = 32;    // Number of explored thetas
	const double dTheta    = (theta_max - theta0) / nTheta;


	/* +---------------------------------+
	 * | STEP 1: INTEGRATING UP TO x = 1 |
	 * +---------------------------------+ */

	const double x_max  = 1.0;
	const double x_over = 0.1;

	// ------------ SHOOTING ------------
	out.open("data/shooting.csv");
	if (!out) exit(5);

	out << "t,x,y,u,v,theta" << endl;  // Output file header
	gTheta = theta0;
	for (int j = 0; j < nTheta; j++) {
		t          = t0;
		double y[] = {0.0, 0.0, v0 * cos(gTheta),
		              v0 * sin(gTheta)};  // Initialize solution
		out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
			<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
			<< endl;  // Output initial condition
		// Integration
		int counter = 0;
		while (y[0] < x_max + x_over && counter < nStep) {
			pVerlet(t, y, RHS, dt, nEq);
			t += dt;

			out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
				<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
				<< endl;
			counter++;
		}
		gTheta += dTheta;
	}
	out.close();


	/* +-----------------------------------+
	 * | STEP 2: FINDING ROOTS OF x(t) = 1 |		Maybe step 1 can be implemented inside x_t
	 * +-----------------------------------+ */

	const double t_low = 1.0;     //!< Lower bound for root searching
	const double t_upp = 2.0;     //!< Upper bound for root searching
	const double eps_t = 1.0e-7;  //!< Precision for root searching

	out.open("data/x_t.csv");
	if (!out) exit(5);

	out << "t*,theta" << endl;
	gTheta = theta0;
	double
		t_star[64];  //!< Array with the times at which x == 1. Length == nTheta
	for (int i = 0; i < nTheta; i++) {
		double roots[4];
		int nRoots;
		try {
			findRoots(x_t, t_low, t_upp, eps_t, roots, nRoots);
		}
		catch (std::runtime_error e) {
			cerr << "You fucked up (step 2). " << e.what() << endl;
		}

		t_star[i] = roots[0];

		out << tau * t_star[i] << "," << gTheta << endl;

		gTheta += dTheta;
	}
	out.close();

	// cout << "t* = ";
	// printVector(t_star, nTheta);


	/* +---------------------------------+
	 * | STEP 3: RESIDUAL PLOT OF y(x=1) |
	 * +---------------------------------+ */

	// ------------ SHOOTING ------------
	std::ofstream out2;
	out2.open("data/shooting2.csv");
	if (!out2) exit(5);

	out2 << "t,x,y,u,v,theta" << endl;
	out2.close();

	out.open("data/residual.csv");
	if (!out) exit(5);

	out << "t*,y*,theta" << endl;  // Output file header
	gTheta = theta0;
	for (int j = 0; j < nTheta; j++) {
		out << tau * t_star[j] << "," << chi * y_t(t_star[j]) << "," << gTheta
			<< endl;
		gTheta += dTheta;
	}
	out.close();


	/* +---------------------------------+
	 * | STEP 4: FINDING ROOTS OF y(x=1) |
	 * +---------------------------------+ */

	// const double theta_low = theta0;     //!< Lower bound for root searching
	// const double theta_upp = theta_max;  //!< Upper bound for root searching
	// const double eps_theta = 1.0e-7;     //!< Precision for root searching

	// out.open("data/optimal.csv");
	// if (!out) exit(5);

	// out << "theta*" << endl;
	// double theta_star[64];  //!< Array with the thetas at which y(x=1) = 0
	// int nTheta_star;        //!< Length of theta_star
	// try {
	// 	// This can't work: y_t is a function of time alone, and I'm looking for
	// 	// roots in theta.
	// 	findRoots(y_t, theta_low, theta_upp, eps_theta, theta_star, nTheta_star, 2);
	// }
	// catch (std::runtime_error e) {
	// 	cerr << "You fucked up (step 4). " << e.what() << endl;
	// }
	// for (int i = 0; i < nTheta_star; i++) {
	// 	out << theta_star[i] << endl;
	// }
	// out.close();

	// cout << "Optimal thetas: ";
	// printVector(theta_star, nTheta_star);


	/* +---------------------------------+
	 * | STEP 5: INTEGRATING WITH THETA* |
	 * +---------------------------------+ */

	out.open("data/final_trajectories.csv");
	if (!out) exit(5);

	out << "t,x,y,u,v,theta" << endl;  // Output file header
	gTheta     = 0.68905;
	t          = t0;
	double y[] = {0.0, 0.0, v0 * cos(gTheta),
					v0 * sin(gTheta)};  // Initialize solution
	out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
		<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
		<< endl;  // Output initial condition
	// Integration
	int counter = 0;
	while (y[0] < x_max && counter < nStep) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
			<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
			<< endl;
		counter++;
	}
	out.close();

	return 0;
}

void RHS(double Y[], double R[]) {
	// double x = Y[0];
	// double y = Y[1];
	double u = Y[2];
	double v = Y[3];

	double mod_v = sqrt(u * u + v * v);

	R[0] = u;
	R[1] = v;
	R[2] = -b * u * mod_v;
	R[3] = -1.0 - b * u * mod_v;
}

double x_t(const double& time) {
	// Integration range
	const double t0    = 0.0;
	const double t_max = time;
	const int nStep    = 1000;
	const double dt    = (t_max - t0) / nStep;

	// Initialisation
	double t      = t0;
	double y[]    = {0.0, 0.0, v0 * cos(gTheta), v0 * sin(gTheta)};
	const int nEq =
		static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	// Integration
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;
	}
	return y[0] - 1.0;
}

double y_t(const double& time) {
	std::ofstream out;
	out.open("data/shooting2.csv", std::ios_base::app);
	if (!out) exit(5);

	// Integration range
	const double t0    = 0.0;
	const double t_max = time;
	const int nStep    = 1000;
	const double dt    = (t_max - t0) / nStep;

	// Initialisation
	double t   = t0;
	double y[] = {0.0, 0.0, v0 * cos(gTheta), v0 * sin(gTheta)};
	const int nEq =
		static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	// Integration
	out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
		<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
		<< endl;  // Output initial condition
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
			<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
			<< endl;
	}
	return y[1];
}
