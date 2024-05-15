/**
 * @file main.cpp
 *
 * @author Francesco Marchisotti
 *
 * @brief Main file for the course's final project.
 *
 * @date 2024-05-15
 */

#include "../../Libs/include/exception.hpp"
#include "../../Libs/include/lin_alg.hpp"
#include "../../Libs/include/ode_solver.hpp"
#include "../../Libs/include/root_finder.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

int numIntegrations = 0;  //!< Number of integrations of the ODEs performed

// Problem data
const static double B      = 4.0e-5;  //!< Drag coefficient [kg/m]
const static double V0     = 9.90;    //!< Initial velocity [m/s]
const static double L      = 10.0;    //!< Target distance  [m]
const static double Y_targ = -0.2;    //!< Target height    [m]
double gTheta;                        //!< Initial launch angle

// Dimensional factors
const static double chi = L;              //!< Space dimensional factor [m]
const static double mu  = 1.0;            //!< Mass dimensional factor [kg]
const static double g   = 9.81;           //!< Gravity [m/s^2]
const static double tau = sqrt(chi / g);  //!< Time dimensional factor [s]

const static double b      = B * chi / mu;    //!< Adimensional friction
const static double v0     = V0 * tau / chi;  //!< Adimensional speed
const static double y_targ = Y_targ / L;      //!< Adimensional target height

/**
 * @brief      Right Hand Side of the system of ODEs.
 *
 * @param[in]  Y  The input values of the variables.
 * @param[out] R  The output values of the variables.
 */
void RHS(const double& t, double Y[], double R[]);

/**
 * @brief     This function returns the linear interpolation between two points
 *            evaluated at a certain x.
 *
 * @param[in] x   The point at which to evaluate the interpolation.
 * @param[in] x1  The x-coordinate of the first point.
 * @param[in] y1  The y-coordinate of the first point.
 * @param[in] x2  The x-coordinate of the second point.
 * @param[in] y2  The y-coordinate of the second point.
 *
 * @return    The interpolated line evaluated at x.
 */
double linearInterp(const double& x, const double& x1, const double& y1,
                    const double& x2, const double& y2);

/**
 * @brief          Function that performs the integration and prints to file.
 *
 * @param[in]      RHSFunc     Right Hand Side of the system of ODEs.
 * @param[in, out] y           Array with the variables. Should be already
 *                             initialised.
 * @param[in]      nEq         Number of equations in the system.
 * @param[in]      t           Starting time of the integration.
 * @param[in]      dt          Time step size.
 * @param[in]      theta       Launch angle.
 * @param[out]     xLast       x-position at penultimate step in time.
 * @param[out]     yLast       y-position at penultimate step in time.
 * @param[in]      maxStep     Maximum number of integration steps.
 * @param[in]      outFile  Output file.
 */
void integration(void (*RHSFunc)(const double& t, double* Y, double* RHS),
                 double y[], const int& nEq, double t, const double& dt,
                 const double& theta, double& xLast, double& yLast,
                 const int& maxStep, std::ofstream& outFile);

/**
 * @brief Prints problem data and dimensional constants to file.
 */
void printConstants();

/**
 * @brief     Generates data for shooting plot.
 *
 * @param[in] thetaMin  Lower bound for launch angle.
 * @param[in] thetaMax  Upper bound for launch angle.
 * @param[in] nTheta    Number of launch angles explored.
 */
void shootingPlot(const double& thetaMin, const double& thetaMax,
                  const double& nTheta);

/**
 * @brief          Integration interface.
 *
 * @param[in,out]  y        Array with the variables.
 * @param[in]      theta    Launch angle.
 * @param[out]     xLast    x-position at penultimate step in time.
 * @param[out]     yLast    y-position at penultimate step in time.
 * @param[in]      outFile  Output file.
 */
void integrate(double y[], const double& theta, double& xLast, double& yLast,
               std::ofstream& outFile);

/**
 * @overload
 *
 * @brief          Integration interface (without file output).
 *
 * @param[in,out]  y      Array with the variables.
 * @param[in]      theta  Launch angle.
 * @param[out]     xLast  x-position at penultimate step in time.
 * @param[out]     yLast  y-position at penultimate step in time.
 */
void integrate(double y[], const double& theta, double& xLast, double& yLast);

/**
 * @overload
 *
 * @brief          Integration interface (without penultimate values).
 *
 * @param[in,out]  y      Array with the variables.
 * @param[in]      theta  Launch angle.
 * @param[in]      outFile  Output file.
 */
void integrate(double y[], const double& theta, std::ofstream& outFile);

/**
 * @overload
 *
 * @brief          Integration interface (without penultimate values and file
 * 				   output).
 *
 * @param[in,out]  y      Array with the variables.
 * @param[in]      theta  Launch angle.
 */
void integrate(double y[], const double& theta);

/**
 * @brief     Residual function for the BVP.
 *
 * @param[in] theta  Launch angle.
 *
 * @return    y(x = 1) - y_targ
 */
double Residual(const double& theta);

int main() {
	printConstants();

	const double thetaMin = 0.65;    // Minimum launch angle
	const double thetaMax = 0.92;    // Maximum launch angle
	const int nTheta      = 32;      // Number of launch angles explored
	const double thetaTol = 1.0e-7;  // Tolerance for root searching

	shootingPlot(thetaMin, thetaMax, nTheta);

	numIntegrations = 0;
	double roots[4];
	int nRoots = -1;
	try {
		findRoots(Residual, thetaMin, thetaMax, thetaTol, roots, nRoots, 2,
		          "secant");

		cout << "Integrations performed: " << numIntegrations << endl;

		cout << "Optimal thetas [rad]: ";
		printVector(roots, nRoots);
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not recognise the error." << endl;
	}

	std::ofstream finalTrajectories;
	finalTrajectories.open("data/finalTrajectories.csv");
	// finalTrajectories.open("data/noFriction.csv");
	finalTrajectories << "t,x,y,u,v,theta" << endl;

	for (int i = 0; i < nRoots; i++) {
		double y[4];
		integrate(y, roots[i], finalTrajectories);
	}

	finalTrajectories.close();
	return 0;
}

double linearInterp(const double& x, const double& x1, const double& y1,
                    const double& x2, const double& y2) {
	double m = (y2 - y1) / (x2 - x1);
	double q = y1 - m * x1;
	double y = m * x + q;
	return y;
}

void printConstants() {
	std::ofstream out;
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
		cerr << "Sorry, could not recognise the error." << endl;
	}

	out.close();
}

void shootingPlot(const double& thetaMin, const double& thetaMax,
                  const double& nTheta) {
	const double dTheta = (thetaMax - thetaMin) / nTheta;

	std::ofstream shooting, residual;
	shooting.open("data/shooting.csv");
	residual.open("data/residual.csv");
	shooting << "t,x,y,u,v,theta" << endl;  // Output csv header
	residual << "theta,res" << endl;        // Output csv header
	for (int i = 0; i < nTheta; i++) {
		double theta = thetaMin + i * dTheta;
		double y[4];

		integrate(y, theta, shooting);

		residual << theta << "," << Residual(theta) << endl;
	}
	shooting.close();
	residual.close();
}

void integrate(double y[], const double& theta, double& xLast, double& yLast,
               std::ofstream& outFile) {
	const double y0[] = {0.0, 0.0, v0 * cos(theta), v0 * sin(theta)};
	const int nEq =
		static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));
	for (int i = 0; i < nEq; i++) y[i] = y0[i];

	double t0         = 0.0;
	const double dt   = 1.0e-2;
	const int maxStep = 10000;

	integration(RHS, y, nEq, t0, dt, theta, xLast, yLast, maxStep, outFile);
}

void integrate(double y[], const double& theta, double& xLast, double& yLast) {
	std::ofstream dummyOutfile;
	integrate(y, theta, xLast, yLast, dummyOutfile);
}

void integrate(double y[], const double& theta, std::ofstream& outFile) {
	double dummyLast;
	integrate(y, theta, dummyLast, dummyLast, outFile);
}

void integrate(double y[], const double& theta) {
	std::ofstream dummyOutfile;
	double dummyLast;
	integrate(y, theta, dummyLast, dummyLast, dummyOutfile);
}

void integration(void (*RHSFunc)(const double& t, double* Y, double* RHS),
                 double y[], const int& nEq, double t, const double& dt,
                 const double& theta, double& xLast, double& yLast,
                 const int& maxStep, std::ofstream& outFile) {
	xLast = y[0];  // Store x-position at previous time step
	yLast = y[1];  // Store y-position at previous time step

	outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
			<< "," << theta << endl;


	numIntegrations++;
	int stepCounter    = 0;
	bool exitCondition = false;
	while (stepCounter < maxStep && !exitCondition) {
		xLast = y[0];
		yLast = y[1];

		rk4Step(t, y, RHSFunc, dt, nEq);
		t += dt;
		stepCounter++;

		outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
				<< "," << theta << endl;

		if (xLast < 1.0 && y[0] > 1.0) exitCondition = true;
	}
}

void RHS(const double& t, double Y[], double R[]) {
	double u = Y[2];
	double v = Y[3];

	double mod_v = sqrt(u * u + v * v);

	R[0] = u;
	R[1] = v;
	R[2] = -b * u * mod_v;
	R[3] = -1.0 - b * u * mod_v;
}

double Residual(const double& theta) {
	double y[4];
	double xLast, yLast;
	integrate(y, theta, xLast, yLast);

	double xCurrent = y[0], yCurrent = y[1];

	return linearInterp(1.0, xLast, yLast, xCurrent, yCurrent) - y_targ;
}
