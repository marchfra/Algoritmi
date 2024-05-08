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
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder_param.hpp"

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
const static double V0 = 1.01;    //!< Initial velocity
const static double L  = 1.0;    //!< Target distance
double gTheta;                    //!< Initial launch angle

// Dimensional factors
const static double chi = L;
const static double mu  = 1.0;
const static double g   = 1.0;
const static double tau = sqrt(chi / g);
const static double b   = B * chi / mu;
const static double v0  = V0 * tau / chi;

/**
 * @brief Right Hand Side of the ODEs system
 *
 * @param[in] Y The input values of the variables.
 * @param[out] R The output values of the variables.
 */
void RHS(double Y[], double R[]);

/**
 * @brief Residual function of x(t)
 *
 * This function takes in input a time value and integrates the system of ODEs
 * upto (and including) that time, and returns the value of x(t) - 1 in order to
 * find the time at which x == 1.
 *
 * @param[in] time The time at which to evaluate x.
 *
 * @return x(t) - 1
 */
double x_t_plot(const double& time);

/**
 * @overload
 *
 * @brief Residual function of x(t, theta)
 *
 * This function takes in input a time value and an initial launch angle and integrates the system of ODEs
 * upto (and including) that time, and returns the value of x(t) - 1 in order to
 * find the time at which x == 1.
 *
 * @param[in] time  The time at which to evaluate x.
 * @param[in] theta The initial launch angle.
 *
 * @return x(t, theta) - 1
 */
double x_t(const double& time, const double& theta);

/**
 * @brief Residual function of y(t)
 *
 * This function takes in input a time value and integrates the system of ODEs
 * upto (and including) that time, and returns the value of y(t*) in order to
 * find the x-position at which y == 0.
 *
 * @param[in] time The time at which to evaluate y.
 *
 * @return y(t)
 */
double y_t_plot(const double& time);

/**
 * @brief Residual function of y(t)
 *
 * This function takes in input a time value and integrates the system of ODEs
 * upto (and including) that time, and returns the value of y(t*) in order to
 * find the x-position at which y == 0.
 *
 * @param[in] time The time at which to evaluate y.
 *
 * @return y(t)
 */
double y_t(const double& time, const double& theta);

/**
 * @brief Residual function.
 *
 * Find the roots of this function to estimate optimal launch angle.
 *
 * @param[in] theta The angle at which to launch the projectile.
 *
 * @return The height of the projectile at x = L
 */
double Residual(const double& theta);

double parabola(double x);

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
	const double theta0    = 0.65;  // Initial guess for launch angle
	const double theta_max = 0.92;  // Max value for theta
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
		// if (j == 0) cout << "!DEBUG  --  " << "y(x + eps) = " << y[1] << endl;
		gTheta += dTheta;
	}
	out.close();


	// -------------------------------------------------------------------------
	// out.open("data/testing.csv");
	// out << "t,x,y,u,v,theta,param" << endl;
	// out.close();
	// try {
	// 	cout << "Expected value       =  1.24443\n";
	// 	gTheta = theta0;
	// 	double ROOTS[] = {-1.0, -2.0, -3.0, -4.0};
	// 	int NROOTS = 1;
	// 	findRoots(x_t_plot, 1.0, 2.0, 1.0e-7, ROOTS, NROOTS);
	// 	// bisection(x_t_plot, 1.0, 2.0, 1.0e-7, ROOTS[0]);
	// 	// cout << "x(t = 1.5) = " << x_t_plot(1.5) << endl;
	// 	cout << "Roots of x(t)        = ";
	// 	printVector(ROOTS, NROOTS);
	// 	findRoots(x_t, gTheta, 1.0, 2.0, 1.0e-7, ROOTS, NROOTS);
	// 	// int NTRY;
	// 	// bisection(x_t, gTheta, 1.0, 2.0, 1.0e-7, 1.0e-7, ROOTS[0], NTRY);
	// 	// cout << "x(t = 1.5, theta) = " << x_t(1.5, gTheta) << endl;
	// 	cout << "Roots of x(t, theta) = ";
	// 	printVector(ROOTS, NROOTS);
	// }
	// catch (...) {
	// 	cout << "Exception" << endl;
	// }
	// -------------------------------------------------------------------------


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
	if (nTheta > 64) throw "nTheta should be at most 64. Lower nTheta or increase size of t_star";
	for (int i = 0; i < nTheta; i++) {
		double roots[4];
		int nRoots;
		try {
			findRoots(x_t_plot, t_low, t_upp, eps_t, roots, nRoots);
		}
		catch (std::runtime_error& e) {
			cerr << "You fucked up (step 2). " << e.what() << endl;
		}

		t_star[i] = roots[0];

		out << tau * t_star[i] << "," << gTheta << endl;

		gTheta += dTheta;
	}
	out.close();

	// cout << "!DEBUG  --  ";
	// cout << "t* = ";
	// printVector(t_star, nTheta);

	// cout << "!DEBUG  --  " << "correct t*(theta = " << theta0 << ") = " << t_star[0] << endl << endl;

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
		double Y = y_t_plot(t_star[j]);
		out << tau * t_star[j] << "," << chi * Y << "," << gTheta
			<< endl;
		cout << "(theta, y*) = (" << gTheta << "," << Y << ")" << endl;
		gTheta += dTheta;
	}
	out.close();


	/* +---------------------------------+
	 * | STEP 4: FINDING ROOTS OF y(x=1) |
	 * +---------------------------------+ */

	const double theta_low = theta0;     //!< Lower bound for root searching
	const double theta_upp = theta_max;  //!< Upper bound for root searching
	const double eps_theta = 1.0e-7;     //!< Precision for root searching

	out.open("data/optimal.csv");
	if (!out) exit(5);

	out << "theta*" << endl;
	double theta_star[4];  //!< Array with the thetas at which y(x=1) = 0
	int nTheta_star = 1;       //!< Length of theta_star
	// double THETA;

	// cout << endl;
	// Residual(theta_low);  // !DEBUG
	// Residual(theta_upp - dTheta);  // !DEBUG
	// cout << endl;

	// cout << endl;
	// double THETAS[] = {theta_low, theta_upp - dTheta};
	// double T_STARS[] = {1.24374, 1.6165};
	// double EXPECTEDS_T[] = {t_star[0], t_star[nTheta - 1]};
	// double EXPECTEDS_Y[] = {-0.0132444, -0.0159998};
	// for (int i = 0; i < 2; i++) {
	// 	cout << "theta    = " << THETAS[i] << endl;
	// 	cout << "t*       = " << T_STARS[i] << endl;
	// 	cout << "Expected = " << EXPECTEDS_T[i] << endl;
	// 	cout << "y(t*)    = " << y_t(T_STARS[i], THETAS[i]) << endl;
	// 	cout << "Expected = " << EXPECTEDS_Y[i] << endl << endl;
	// }

	try {
		std::ofstream out2;
		out2.open("data/optimal_search.csv");
		out2 << "t*,y*,theta" << endl;
		out2.close();

		findRoots(Residual, theta_low, theta_upp, eps_theta, theta_star, nTheta_star);
		// bisection(Residual, theta_low, theta_upp, eps_theta, THETA);
		// cout << "\n\n\n The root is in theta = " << THETA << "\n\n";
		if (nTheta_star > 2) throw "There should be a maximum of 2 solutions.";
	}
	catch (std::runtime_error e) {
		cerr << "You fucked up (step 4). " << e.what() << endl;
	}
	catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	}
	catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}
	// for (int i = 0; i < nTheta_star; i++) theta_star[i] = THETA;
	for (int i = 0; i < nTheta_star; i++) {
		out << theta_star[i] << endl;
	}
	out.close();

	cout << "Optimal thetas: ";
	printVector(theta_star, nTheta_star);


	/* +---------------------------------+
	 * | STEP 5: INTEGRATING WITH THETA* |
	 * +---------------------------------+ */

	out.open("data/final_trajectories.csv");
	if (!out) exit(5);

	out << "t,x,y,u,v,theta" << endl;  // Output file header
	for (int i = 0; i < nTheta_star; i++) {
		t          = t0;
		gTheta     = theta_star[i];
		double y[] = {0.0, 0.0, v0 * cos(gTheta),
						v0 * sin(gTheta)};  // Initialize solution
		out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
			<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
			<< endl;  // Output initial condition
		// Integration
		int ctr = 0;
		while (y[0] < x_max && ctr < nStep) {
			pVerlet(t, y, RHS, dt, nEq);
			t += dt;

			out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
				<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
				<< endl;
			ctr++;
		}
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

double x_t_plot(const double& time) {
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

	// std::ofstream out;
	// out.open("data/testing.csv", std::ios_base::app);

	// out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << gTheta << ",0" << endl;
	// Integration
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		// out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << gTheta << ",0" << endl;
	}

	// out.close();
	return y[0] - 1.0;
}

double x_t(const double& time, const double& theta) {
	// Integration range
	const double t0    = 0.0;
	const double t_max = time;
	const int nStep    = 1000;
	const double dt    = (t_max - t0) / nStep;

	// Initialisation
	double t      = t0;
	double y[]    = {0.0, 0.0, v0 * cos(theta), v0 * sin(theta)};
	const int nEq =
		static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	// std::ofstream out;
	// out.open("data/testing.csv", std::ios_base::app);

	// out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta << ",1" << endl;
	// Integration
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		// out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << theta << ",1" << endl;
	}

	// out.close();
	return y[0] - 1.0;
}

double y_t_plot(const double& time) {
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
	out.close();
	return y[1];
}

double y_t(const double& time, const double& theta) {
	// Integration range
	const double t0    = 0.0;
	const double t_max = time;
	const int nStep    = 1000;
	const double dt    = (t_max - t0) / nStep;

	// Initialisation
	double t   = t0;
	double y[] = {0.0, 0.0, v0 * cos(theta), v0 * sin(theta)};
	const int nEq =
		static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));

	// Integration
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;
	}
	return y[1];
}

double Residual(const double& theta) {
	cout << "!DEBUG  --  " << "theta = " << theta << endl;
	std::ofstream out;

	// Integration range
	const double t0    = 0.0;                // Starting point
	const double t_max = 10.0;               // Ending point
	const double dt    = 0.01;               // Time increment
	const int nStep    = (t_max - t0) / dt;  // Number of steps

	// Initialisation
	double t          = t0;
	const double y0[] = {0.0, 0.0, v0 * cos(theta), v0 * sin(theta)};
	const int nEq =
		static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));
	double y[nEq];


	/* +---------------------------------+
	 * | STEP 1: INTEGRATING UP TO x = 1 |
	 * +---------------------------------+ */

	const double x_max  = 1.0;
	const double x_over = 0.1;

	t          = t0;
	for (int i = 0; i < nEq; i++) y[i] = y0[i];  // Initialize solution
	// Integration
	int ctr = 0;
	while (y[0] < x_max + x_over && ctr < nStep) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		ctr++;
	}
	cout << "!DEBUG  --  " << "y(x + eps) = " << y[1] << endl;


	/* +-----------------------------------+
	 * | STEP 2: FINDING ROOTS OF x(t) = 1 |		Maybe step 1 can be implemented inside x_t
	 * +-----------------------------------+ */

	const double t_low = 1.0;     //!< Lower bound for root searching
	const double t_upp = 2.0;     //!< Upper bound for root searching
	const double eps_t = 1.0e-7;  //!< Precision for root searching

	double roots[4];
	int nRoots;
	try {
		findRoots(x_t, theta, t_low, t_upp, eps_t, roots, nRoots);

		// cout << "!DEBUG  --  " << "Root of x(t) - 1 for theta = " << theta << " is t* = " << roots[0] << endl;

		if (nRoots == 0) throw "x doesn't cross the target. Increase initial angle";
		else if (nRoots > 1) throw "x crosses the target at two points in time";
	}
	catch (std::runtime_error e) {
		cerr << "You fucked up (step 2). " << e.what() << endl;
	}
	catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	}
	catch (...) {
		cerr << "Sorry, could not catch the error, whatever it is." << endl;
	}

	const double t_star = roots[0];
	cout << "!DEBUG  --  " << "t* = " << tau * t_star << endl;

	double y_star = y_t(t_star, theta);
	cout << "!DEBUG  --  " << "y* = " << chi * y_star << endl << endl;

	out.open("data/optimal_search.csv", std::ios_base::app);
	out << tau * t_star << "," << chi * y_star/*parabola(theta)*/ << "," << theta << endl;
	out.close();

	return y_star;
}

double parabola(double x) {
	return -20.0 * (x - 0.787)*(x - 0.787) + 0.19;
}
