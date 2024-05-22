/**
 * @file main.cpp
 *
 * @author Francesco Marchisotti
 *
 * @brief Main file for the course's final project.
 *
 * @date 2024-05-15
 *
 * @todo quadraticInterp using Gaussian elimination to solve the system of
 * equations
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

#define FRICTION 0
const static int gOrder = 1;

int numIntegrations = 0;  //!< Number of integrations of the ODEs performed

double g_dt = 1.0e-3;

// Problem data
#if FRICTION
const static double B      = 4.0e-5;  //!< Drag coefficient [kg/m]
const static double V0     = 9.90;    //!< Initial velocity [m/s]
const static double L      = 10.0;    //!< Target distance  [m]
const static double Y_targ = -0.2;    //!< Target height    [m]
double gTheta;                        //!< Initial launch angle
#else
const static double B      = 0.0;   //!< Drag coefficient [kg/m]
const static double V0     = 10.0;  //!< Initial velocity [m/s]
const static double L      = 10.0;  //!< Target distance  [m]
const static double Y_targ = 0.0;   //!< Target height    [m]
double gTheta;                      //!< Initial launch angle
#endif

// Dimensional factors
const static double chi = L;              //!< Space dimensional factor [m]
const static double mu  = 1.0;            //!< Mass dimensional factor [kg]
const static double g   = 9.81;           //!< Gravity [m/s^2]
const static double tau = sqrt(chi / g);  //!< Time dimensional factor [s]

const static double b      = B * chi / mu;    //!< Adimensional friction
const static double v0     = V0 * tau / chi;  //!< Adimensional speed
const static double x_targ = 1.0;             //!< Adimensional target distance
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
 * @brief     This function returns the quadratic interpolation between three
 * 			  points evaluated at a certain x.
 *
 * @param[in] x   The point at which to evaluate the interpolation.
 * @param[in] x1  The x-coordinate of the first point.
 * @param[in] y1  The y-coordinate of the first point.
 * @param[in] x2  The x-coordinate of the second point.
 * @param[in] y2  The y-coordinate of the second point.
 * @param[in] x3  The x-coordinate of the third point.
 * @param[in] y3  The y-coordinate of the third point.
 *
 * @return    The interpolated parabola evaluated at x.
 */
double quadraticInterp(const double& x, const double& x1, const double& y1,
                       const double& x2, const double& y2, const double& x3,
                       const double& y3);

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
 * @param[out]     xLast       Array with x-position at previous time step.
 * @param[out]     yLast       Array with y-position at previous time step.
 * @param[in]      order       Numer of previous times to save.
 * @param[in]      maxStep     Maximum number of integration steps.
 * @param[in]      outFile  Output file.
 */
void integration(void (*RHSFunc)(const double& t, double* Y, double* RHS),
                 double y[], const int& nEq, double t, const double& dt,
                 const double& theta, double xLast[], double yLast[],
                 const int& order, const int& maxStep, std::ofstream& outFile);

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
 * @param[out]     xLast    Array with x-position at previous time step.
 * @param[out]     yLast    Array with y-position at previous time step.
 * @param[in]      order    Numer of previous times to save.
 * @param[in]      outFile  Output file.
 */
void integrate(double y[], const double& theta, double xLast[], double yLast[],
               const int& order, std::ofstream& outFile);

/**
 * @overload
 *
 * @brief          Integration interface (without file output).
 *
 * @param[in,out]  y      Array with the variables.
 * @param[in]      theta  Launch angle.
 * @param[out]     xLast  Array with x-position at previous time step.
 * @param[out]     yLast  Array with y-position at previous time step.
 * @param[in]      order  Numer of previous times to save.
 */
void integrate(double y[], const double& theta, double& xLast, double& yLast);

/**
 * @overload
 *
 * @brief          Integration interface (without past values).
 *
 * @param[in,out]  y        Array with the variables.
 * @param[in]      theta    Launch angle.
 * @param[in]      outFile  Output file.
 */
void integrate(double y[], const double& theta, std::ofstream& outFile);

/**
 * @overload
 *
 * @brief          Integration interface (without past values and file output).
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
#if FRICTION
	cout << "===== FRICTION =====" << endl << endl;
#else
	cout << "===== NO FRICTION =====" << endl << endl;
#endif

	printConstants();

	const double thetaMin = 0.65;    // Minimum launch angle
	const double thetaMax = 0.92;    // Maximum launch angle
	const int nTheta      = 32;      // Number of launch angles explored
	const double thetaTol = 1.0e-7;  // Tolerance for root searching

	// std::ofstream linInterp;
	// linInterp.open("data/linear.csv");
	// linInterp << "x,y" << endl;
	// for (int i = 0; i < 100; i++) {
	// 	double x = -10.0 + 0.2 * i;
	// 	linInterp << x << "," << linearInterp(x, -3, 4, 2, 1.5) << endl;
	// }
	// linInterp.close();

	// std::ofstream quadrInterp;
	// quadrInterp.open("data/quadratic.csv");
	// quadrInterp << "x,y" << endl;
	// for (int i = 0; i < 100; i++) {
	// 	double x = -10.0 + 0.2 * i;
	// 	quadrInterp << x << "," << quadraticInterp(x, -3, 4, 2, 1.5, 1.1, 2.3) << endl;
	// }
	// quadrInterp.close();

	shootingPlot(thetaMin, thetaMax, nTheta);

	numIntegrations = 0;
	double roots[4];
	int nRoots = -1;
	try {
		findRoots(Residual, thetaMin, thetaMax, thetaTol, roots, nRoots, 4,
		          "secant");

		cout << "Integrations performed: " << numIntegrations << endl;

		cout.precision(7);
		cout << "Optimal thetas [rad] (+/- " << thetaTol << "): ";
		printVector(roots, nRoots);
	} catch (std::exception& err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not recognise the error." << endl;
	}

	/* +-------------+
	 * | Convergence |
	 * +-------------+ */

	// std::ofstream out2;
	// out2.open("data/convergence.csv", std::ios_base::app);
	// out2 << "dt,theta1,theta2" << endl;
	// out2 << g_dt << "," << roots[0] << "," << roots[1] << endl;

	// for (int i = 0; i < nRoots; i++) {
	// 	out2 << g_dt << "," << Residual(roots[i]) << "," << roots[i] << endl;
	// }

	// out2.close();


	std::ofstream finalTrajectories;
#if FRICTION
	finalTrajectories.open("data/finalTrajectories.csv");
#else
	finalTrajectories.open("data/noFriction.csv");
#endif
	finalTrajectories << "t,x,y,u,v,theta" << endl;

	// roots[0] = 0.6877629863480418;
	// roots[1] = 0.8830333404468548;

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

double quadraticInterp(const double& x, const double& x1, const double& y1,
                       const double& x2, const double& y2, const double& x3,
                       const double& y3) {
	const int nPoints = 3;

	double** M;
	M    = new double*[nPoints];
	M[0] = new double[nPoints * nPoints];
	for (int i = 1; i < nPoints; i++) M[i] = M[i - 1] + nPoints;

	// Define coefficient matrix
	M[0][0] = x1 * x1;
	M[0][1] = x1;
	M[0][2] = 1;
	M[1][0] = x2 * x2;
	M[1][1] = x2;
	M[1][2] = 1;
	M[2][0] = x3 * x3;
	M[2][1] = x3;
	M[2][2] = 1;

	double v[nPoints] = {y1, y2, y3};
	double coeffs[nPoints];  //<! Array with the coefficients of the parabola

	solveLinSystem(M, v, coeffs, nPoints);

	delete[] M[0];
	delete[] M;

	return coeffs[0] * x * x + coeffs[1] * x + coeffs[2];
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

void integrate(double y[], const double& theta, double xLast[], double yLast[],
               const int& order, std::ofstream& outFile) {
	const double y0[] = {0.0, 0.0, v0 * cos(theta), v0 * sin(theta)};
	const int nEq =
		static_cast<int>(sizeof(y0)) / static_cast<int>(sizeof(y0[0]));
	for (int i = 0; i < nEq; i++) y[i] = y0[i];

	double t0         = 0.0;
	const double dt   = g_dt;
	const int maxStep = int(2 / dt);

	integration(RHS, y, nEq, t0, dt, theta, xLast, yLast, order, maxStep,
	            outFile);
}

void integrate(double y[], const double& theta, double xLast[], double yLast[],
               const int& order) {
	std::ofstream dummyOutfile;
	integrate(y, theta, xLast, yLast, order, dummyOutfile);
}

void integrate(double y[], const double& theta, std::ofstream& outFile) {
	double dummyLast[2];
	integrate(y, theta, dummyLast, dummyLast, 1, outFile);
}

void integrate(double y[], const double& theta) {
	std::ofstream dummyOutfile;
	double dummyLast[2];
	integrate(y, theta, dummyLast, dummyLast, 1, dummyOutfile);
}

void integration(void (*RHSFunc)(const double& t, double* Y, double* RHS),
                 double y[], const int& nEq, double t, const double& dt,
                 const double& theta, double xLast[], double yLast[],
                 const int& order, const int& maxStep, std::ofstream& outFile) {
	outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
			<< "," << theta << endl;

	numIntegrations++;
	int stepCounter    = 0;
	bool exitCondition = false;
	while (stepCounter < maxStep && !exitCondition) {
		switch (order) {
		case 1:
			xLast[0] = y[0];  // Store x-position at previous time step
			yLast[0] = y[1];  // Store y-position at previous time step
			break;
		case 2:
			xLast[0] = xLast[1];
			yLast[0] = yLast[1];
			xLast[1] = y[0];
			yLast[1] = y[0];
			break;
		default: throw exception("order must be either 1 or 2"); break;
		}

		rk4Step(t, y, RHSFunc, dt, nEq);
		t += dt;
		stepCounter++;

		outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
				<< "," << theta << endl;

		if (xLast[0] < x_targ && y[0] > x_targ) exitCondition = true;
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
	double xLast[2], yLast[2];
	integrate(y, theta, xLast, yLast, gOrder);

	double xCurrent = y[0], yCurrent = y[1];

	switch (gOrder) {
	case 1:
		return linearInterp(x_targ, xLast[0], yLast[0], xCurrent, yCurrent) -
		       y_targ;
		break;
	case 2:
		return quadraticInterp(x_targ, xLast[0], yLast[0], xLast[1], yLast[1],
		                       xCurrent, yCurrent) -
		       y_targ;
		break;
	default: throw exception("gOrder must be either 1 or 2."); break;
	}
}
