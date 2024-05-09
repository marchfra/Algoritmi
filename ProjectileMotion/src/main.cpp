/**
 * @file main.cpp
 *
 * @author Francesco Marchisotti
 *
 * @brief Main file for the course's final project
 *
 * @date 2024-05-04
 */

#include "../../Libs/include/exception.hpp"
#include "../../Libs/include/lin_alg.hpp"
#include "../../Libs/include/ode_solver.hpp"
#include "../../Libs/include/root_finder.hpp"
#include "../../Libs/include/root_finder_param.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

// Problem data
const static double B      = 4.0e-5;  //!< Drag coefficient [kg/m]
const static double V0     = 9.90;    //!< Initial velocity [m/s]
const static double L      = 10.0;    //!< Target distance  [m]
const static double Y_targ = -0.2;    //!< Target height    [m]
double gTheta;                        //!< Initial launch angle

// Dimensional factors
const static double chi    = L;               //!< Space dimensional factor
const static double mu     = 1.0;             //!< Mass dimensional factor
const static double g      = 9.81;            //!< Gravity [m/s^2]
const static double tau    = sqrt(chi / g);   //!< Time dimensional factor
const static double b      = B * chi / mu;    //!< Adimensional friction
const static double v0     = V0 * tau / chi;  //!< Adimensional speed
const static double y_targ = Y_targ / L;      //!< Adimensional target height

/**
 * @brief      Right Hand Side of the system of ODEs.
 *
 * @param[in]  Y  The input values of the variables.
 * @param[out] R  The output values of the variables.
 */
void RHS(double Y[], double R[]);

/**
 * @brief     Computes the x-component of the trajectory.
 *
 * This function takes in input a time value and an initial launch angle,
 * integrates the system of ODEs upto (and including) that time, and returns the
 * value of x(t, theta) - 1 in order to find the time t1 at which x(t1, theta)
 * = 1.
 *
 * @param[in] time   The time at which to evaluate x.
 * @param[in] theta  The initial launch angle.
 *
 * @return    x(t, theta) - 1
 */
double x_t(const double& time, const double& theta);

/**
 * @brief     Residual function of y(t).
 *
 * This function takes in input a time value, integrates the system of ODEs up
 * to (and including) that time, and returns the value of y(t) in order to make
 * the residual plot of y(t1, theta). While integrating the equations, if a true
 * print parameter is passed, it also prints to file the values of the dependent
 * variables, in order to make a plot of the trajectory up to x = 1.
 *
 * @param[in] time   The time at which to evaluate y.
 * @param[in] theta  The initial launch angle.
 * @param[in] print  Whether to print to file.
 *
 * @return    y(t, theta)
 */
double y_t(const double& time, const double& theta, const bool print = false);

/**
 * @brief     Residual function of the problem.
 *
 * This function computes y(x = L).
 * Find the roots of this function to estimate optimal launch angle.
 *
 * @param[in] theta  The angle at which to launch the projectile.
 *
 * @return    The height of the projectile at x = L.
 */
double Residual(const double& theta);

/**
 * @brief     Function that handles all the exploration plots.
 *
 * @param[in] t0         Starting time of the integration.
 * @param[in] dt         Time step size.
 * @param[in] nStep      Number of time steps.
 * @param[in] theta_min  Smallest initial angle.
 * @param[in] dTheta     Angle step size.
 * @param[in] nTheta     Number of angles explored.
 */
void plot(const double& t0, const double& dt, const int& nStep,
          const double& theta_min, const double& dTheta, const int& nTheta);

int main() {
	std::ofstream out;  //!< Object containing output file

	// Print constants to file
	try {
		out.open("data/constants.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "chi,tau,mu,B,b,V0,Y_targ" << endl;
		out << chi << "," << tau << "," << mu << "," << B << "," << b << ","
			<< V0 << "," << Y_targ << endl;

		out.close();
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}


	/* +----------------------+
	 * | STEP 0: MAKING PLOTS |
	 * +----------------------+ */

	// Integration range
	const double t0    = 0.0;                   //!< Starting point
	const double t_max = 10.0;                  //!< Ending point
	const int nStep    = 1000;                  //!< Number of steps
	const double dt    = (t_max - t0) / nStep;  //!< Time step size

	// Initial angle
	const double theta_min = 0.65;  //!< Min guess for launch angle
	const double theta_max = 0.92;  //!< Max value for theta
	const int nTheta       = 32;    //!< Number of explored thetas
	const double dTheta =
		(theta_max - theta_min) / nTheta;  //!< Theta step size

	plot(t0, dt, nStep, theta_min, dTheta, nTheta);


	/* +------------------------------------------+
	 * | STEP 1: FINDING ROOTS OF y(x=1) - y_targ |
	 * +------------------------------------------+ */

	const double eps_theta = 1.0e-7;  //!< Precision for root searching

	double theta_star[4];  //!< Array with the thetas at which y(x=1) = 0
	int nTheta_star = 1;   //!< Length of theta_star

	try {
		out.open("data/optimal.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "theta*" << endl;

		std::ofstream out2;
		out2.open("data/optimal_search.csv");
		if (!out2.good()) {
			out.close();
			throw exception("Invalid file.");
		}
		out2 << "t1,y1,theta" << endl;
		out2.close();

		// Find the roots of the residual function
		findRoots(Residual, theta_min, theta_max, eps_theta, theta_star,
		          nTheta_star, 2, "secant");
		for (int i = 0; i < nTheta_star; i++) out << theta_star[i] << endl;
		out.close();

		// Print optimal launch angles
		cout << "Optimal thetas [rad]: ";
		printVector(theta_star, nTheta_star);
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}


	/* +---------------------------------+
	 * | STEP 2: INTEGRATING WITH THETA* |
	 * +---------------------------------+ */

	try {
		out.open("data/final_trajectories.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "t,x,y,u,v,theta" << endl;  // Output file header
		for (int i = 0; i < nTheta_star; i++) {
			double t     = t0;             // Reinitialise time
			double theta = theta_star[i];  // Select one of the optimal angles
			double y[]   = {0.0, 0.0, v0 * cos(theta),
			                v0 * sin(theta)};  // Initialize solution
			int nEq =
				static_cast<int>(sizeof(y)) / static_cast<int>(sizeof(y[0]));
			out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
				<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << theta
				<< endl;  // Output initial condition
			// Integration
			int ctr = 0;  // Eternal loop exit condition
			while (y[0] < 1.0 && ctr < nStep) {
				pVerlet(t, y, RHS, dt, nEq);
				t += dt;

				out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
					<< chi / tau * y[2] << "," << chi / tau * y[3] << ","
					<< theta << endl;
				ctr++;
			}
		}
		out.close();
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}

	return 0;
}

void RHS(double Y[], double R[]) {
	double u = Y[2];
	double v = Y[3];

	double mod_v = sqrt(u * u + v * v);

	R[0] = u;
	R[1] = v;
	R[2] = -b * u * mod_v;
	R[3] = -1.0 - b * u * mod_v;
}

double x_t(const double& time, const double& theta) {
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

	return y[0] - 1.0;
}

double y_t(const double& time, const double& theta, const bool print) {
	std::ofstream out;
	out.open("data/shooting2.csv", std::ios_base::app);
	if (!out.good() && print) throw exception("Invalid file");

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
	if (print)
		out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
			<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << theta
			<< endl;  // Output initial condition
	for (int i = 0; i < nStep; i++) {
		pVerlet(t, y, RHS, dt, nEq);
		t += dt;

		if (print)
			out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
				<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << theta
				<< endl;
	}
	out.close();
	return y[1] - y_targ;
}

double Residual(const double& theta) {
	const double t_low = 1.0;     //!< Lower bound for root searching
	const double t_upp = 2.0;     //!< Upper bound for root searching
	const double eps_t = 1.0e-7;  //!< Precision for root searching

	// Finding roots of x(t) - 1
	double roots[4];
	int nRoots;
	findRoots(x_t, theta, t_low, t_upp, eps_t, roots, nRoots, 2, "secant");

	if (nRoots == 0)
		throw exception("x doesn't cross the target. Increase initial angle.");
	else if (nRoots > 1)
		throw exception("x crosses the target at two points in time.");

	std::ofstream out;
	out.open("data/optimal_search.csv", std::ios_base::app);
	if (!out.good()) throw exception("Invalid file");

	const double t1     = roots[0];        //!< Root of x(t) - 1
	const double y_star = y_t(t1, theta);  //!< Height at t1

	out << tau * t1 << "," << chi * y_star << "," << theta << endl;
	out.close();

	return y_star;
}

void plot(const double& t0, const double& dt, const int& nStep,
          const double& theta_min, const double& dTheta, const int& nTheta) {
	// Initialisation
	double t          = t0;  //!< Time variable
	const double y0[] = {0.0, 0.0, v0 * cos(gTheta),
	                     v0 * sin(gTheta)};  //!< Initialised variables array
	const int nEq     = static_cast<int>(sizeof(y0)) /
	                static_cast<int>(sizeof(y0[0]));  //!< Number of equations


	/* +------------------------------+
	 * | STEP 1: SHOOTING UP TO x = 1 |
	 * +------------------------------+ */

	const double x_max  = 1.0;  //!< x-target
	const double x_over = 0.1;  //!< Overshoot on x-target

	std::ofstream out;

	try {
		out.open("data/shooting.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "t,x,y,u,v,theta" << endl;  // Output file header
		gTheta = theta_min;                // Set theta to min guess
		for (int j = 0; j < nTheta; j++) {
			t          = t0;  // Reinitialise time
			double y[] = {0.0, 0.0, v0 * cos(gTheta),
			              v0 * sin(gTheta)};  // Initialize solution
			out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
				<< chi / tau * y[2] << "," << chi / tau * y[3] << "," << gTheta
				<< endl;  // Output initial condition
			// Integration
			int ctr = 0;  //!< Eternal loop exit exit condition
			while (y[0] < x_max + x_over && ctr < nStep) {
				pVerlet(t, y, RHS, dt, nEq);
				t += dt;

				out << tau * t << "," << chi * y[0] << "," << chi * y[1] << ","
					<< chi / tau * y[2] << "," << chi / tau * y[3] << ","
					<< gTheta << endl;
				ctr++;
			}
			gTheta += dTheta;  // Increment theta
		}
		out.close();
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}


	/* +-----------------------------------+
	 * | STEP 2: FINDING ROOTS OF x(t) - 1 |
	 * +-----------------------------------+ */

	const double t_low = 1.0;     //!< Lower bound for root searching
	const double t_upp = 2.0;     //!< Upper bound for root searching
	const double eps_t = 1.0e-7;  //!< Precision for root searching

	double t1[64];  //!< Array with the times at which x = 1. Length = nTheta
	try {
		if (nTheta > 64)
			throw exception("nTheta should be at most 64. Lower nTheta or "
			                "increase size of t1.");

		out.open("data/x_t.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "t1,theta" << endl;
		gTheta = theta_min;  // Initialise theta
		for (int i = 0; i < nTheta; i++) {
			double roots[4];  //!< Array with the roots of x(t) - 1
			int nRoots;       //!< Number of roots. Should be 1.
			findRoots(x_t, gTheta, t_low, t_upp, eps_t, roots, nRoots, 2,
			          "secant");
			if (nRoots == 0)
				throw exception(
					"x doesn't cross the target. Increase initial angle.");
			else if (nRoots > 1)
				throw exception("x crosses the target at two points in time.");

			t1[i] = roots[0];  // Save the root found inside the array with all
			                   // roots
			out << tau * t1[i] << "," << gTheta << endl;
			gTheta += dTheta;  // Increment theta
		}
		out.close();
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}

	/* +---------------------------------+
	 * | STEP 3: RESIDUAL PLOT OF y(x=1) |
	 * +---------------------------------+ */

	// ------------ SHOOTING ------------
	std::ofstream out2;
	try {
		out2.open("data/shooting2.csv");
		if (!out2.good()) throw exception("Invalid file.");

		out2 << "t,x,y,u,v,theta" << endl;
		out2.close();

		out.open("data/residual.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "t1,y1,theta" << endl;  // Output file header
		gTheta = theta_min;            // Initialise theta
		for (int j = 0; j < nTheta; j++) {
			out << tau * t1[j] << "," << chi * y_t(t1[j], gTheta, true) << ","
				<< gTheta << endl;
			gTheta += dTheta;  // Increment theta
		}
		out.close();
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not catch the error whatever it is." << endl;
	}
}
