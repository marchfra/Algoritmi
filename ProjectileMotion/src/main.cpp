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

#define FRICTION 1
const static int gOrder = 2;  //<! Selects order of polynomial interpolation

int numIntegrations = 0;  //!< Number of integrations of the ODEs performed

double g_dt = 1.0e-3;

// Problem data
#if FRICTION
const static double B     = 4.0e-5;  //!< Drag coefficient [kg/m]
const static double V0    = 9.90;    //!< Initial velocity [m/s]
const static double L     = 10.0;    //!< Target distance  [m]
const static double YTarg = -0.2;    //!< Target height    [m]
double gTheta;                       //!< Initial launch angle
#else
const static double B     = 0.0;   //!< Drag coefficient [kg/m]
const static double V0    = 10.0;  //!< Initial velocity [m/s]
const static double L     = 10.0;  //!< Target distance  [m]
const static double YTarg = 0.0;   //!< Target height    [m]
double gTheta;                     //!< Initial launch angle
#endif

// Dimensional factors
const static double chi = L;              //!< Space dimensional factor [m]
const static double mu  = 1.0;            //!< Mass dimensional factor [kg]
const static double g   = 9.81;           //!< Gravity [m/s^2]
const static double tau = sqrt(chi / g);  //!< Time dimensional factor [s]

const static double b     = B * chi / mu;    //!< Adimensional friction
const static double v0    = V0 * tau / chi;  //!< Adimensional speed
const static double xTarg = 1.0;             //!< Adimensional target distance
const static double yTarg = YTarg / L;       //!< Adimensional target height

/**
 * @brief Prints problem data and dimensional constants to file.
 */
void printConstants();

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
double linearInterp(const double &x, const double &x1, const double &y1,
                    const double &x2, const double &y2);

/**
 * @brief               This function returns the polinomial interpolation
 * 						between a sufficient number of points evaluated at a
 * 						certain x.
 *
 * @param[in] x         The point at which to evaluate the interpolation.
 * @param[in] xLast     The x-coordinates of the first ``order`` points.
 * @param[in] yLast     The y-coordinates of the first ``order`` points.
 * @param[in] xCurrent  The x-coordinate of the last point.
 * @param[in] yCurrent  The y-coordinate of the last point.
 * @param[in] order     The order or the polynomial.
 *
 * @return    The interpolated line evaluated at x.
 */
double polInterp(const double &x, double xLast[], double yLast[],
                 const double &xCurrent, const double &yCurrent,
                 const int &order);

/**
 * @brief Tests polInterp function with degree 1.
 */
void linInterpTest();

/**
 * @brief Tests polInterp function with degree 2.
 */
void quadraticInterpTest();

/**
 * @brief      Right Hand Side of the system of ODEs.
 *
 * @param[in]  Y  The input values of the variables.
 * @param[out] R  The output values of the variables.
 */
void RHS(const double &t, double Y[], double R[]);

/**
 * @brief     Generates data for shooting plot.
 *
 * @param[in] thetaMin  Lower bound for launch angle.
 * @param[in] thetaMax  Upper bound for launch angle.
 * @param[in] nTheta    Number of launch angles explored.
 */
void shootingPlot(const double &thetaMin, const double &thetaMax,
                  const double &nTheta);

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
void integrate(double y[], const double &theta, double xLast[], double yLast[],
               const int &order, std::ofstream &outFile);

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
void integrate(double y[], const double &theta, double &xLast, double &yLast);

/**
 * @overload
 *
 * @brief          Integration interface (without past values).
 *
 * @param[in,out]  y        Array with the variables.
 * @param[in]      theta    Launch angle.
 * @param[in]      outFile  Output file.
 */
void integrate(double y[], const double &theta, std::ofstream &outFile);

/**
 * @overload
 *
 * @brief          Integration interface (without past values and file output).
 *
 * @param[in,out]  y      Array with the variables.
 * @param[in]      theta  Launch angle.
 */
void integrate(double y[], const double &theta);

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
void integration(void (*RHSFunc)(const double &t, double *Y, double *RHS),
                 double y[], const int &nEq, double t, const double &dt,
                 const double &theta, double xLast[], double yLast[],
                 const int &order, const int &maxStep, std::ofstream &outFile);

/**
 * @brief     Residual function for the BVP.
 *
 * @param[in] theta  Launch angle.
 *
 * @return    y(x = 1) - yTarg
 */
double Residual(const double &theta);

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

	linInterpTest();
	quadraticInterpTest();

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
	} catch (std::exception &err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not recognise the error." << endl;
	}

	// std::ofstream interp;
	// interp.open("data/interpolationOrder.csv", std::ios_base::app);
	// for (int i = 0; i < nRoots; i++) {
	// 	interp.precision(17);
	// 	interp << gOrder << "," << roots[i] << "," << i + 1 << endl;
	// }
	// interp.close();


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

void printConstants() {
	std::ofstream out;
	try {
		out.open("data/constants.csv");
		if (!out.good()) throw exception("Invalid file.");

		out << "chi,tau,mu,B,b,V0,Y_targ" << endl;
		out << chi << "," << tau << "," << mu << "," << B << "," << b << ","
			<< V0 << "," << YTarg << endl;

		out.close();
	} catch (std::exception &err) {
		cerr << "Caught " << typeid(err).name() << " : " << err.what() << endl;
	} catch (...) {
		cerr << "Sorry, could not recognise the error." << endl;
	}

	out.close();
}

double linearInterp(const double &x, const double &x1, const double &y1,
                    const double &x2, const double &y2) {
	double m = (y2 - y1) / (x2 - x1);
	double q = y1 - m * x1;
	double y = m * x + q;
	return y;
}

double polInterp(const double &x, double xLast[], double yLast[],
                 const double &xCurrent, const double &yCurrent,
                 const int &order) {
	const int nPoints = order + 1;

	double **M;
	M    = new double *[nPoints];
	M[0] = new double[nPoints * nPoints];
	for (int i = 1; i < nPoints; i++) M[i] = M[i - 1] + nPoints;

	double *v;
	v = new double[nPoints];

	// Define coefficient matrix
	for (int i = 0; i < nPoints; i++) {
		M[i][nPoints - 1] = 1;
		if (i != nPoints - 1) v[i] = yLast[i];
		else v[i] = yCurrent;
		for (int j = nPoints - 2; j >= 0; j--)
			if (i != nPoints - 1) M[i][j] = M[i][j + 1] * xLast[i];
			else M[i][j] = M[i][j + 1] * xCurrent;
	}

	double *coeffs;  //<! Array with the coefficients of the polynomial
	coeffs = new double[nPoints];

	solveLinSystem(M, v, coeffs, nPoints);

	delete[] M[0];
	delete[] M;
	delete[] v;

	double value    = 0.0;
	double powerOfX = 1.0;
	for (int i = nPoints - 1; i >= 0; i--) {
		value += coeffs[i] * powerOfX;
		powerOfX *= x;
	}

	return value;
}

void linInterpTest() {
	std::ofstream linInterp;
	linInterp.open("data/linear.csv");
	linInterp << "x,y" << endl;
	double xL[] = {-3.0}, yL[] = {4.0};
	double xC = 2.0, yC = 1.5;
	for (int i = 0; i < 100; i++) {
		double x = -10.0 + 0.2 * i;
		linInterp << x << "," << polInterp(x, xL, yL, xC, yC, 1) << endl;
	}
	linInterp.close();
}

void quadraticInterpTest() {
	std::ofstream quadrInterp;
	quadrInterp.open("data/quadratic.csv");
	quadrInterp << "x,y" << endl;
	double xL[] = {-3.0, 2.0}, yL[] = {4.0, 1.5};
	double xC = 1.1, yC = 2.3;
	for (int i = 0; i < 100; i++) {
		double x = -10.0 + 0.2 * i;
		quadrInterp << x << "," << polInterp(x, xL, yL, xC, yC, 2) << endl;
	}
	quadrInterp.close();
}

void RHS(const double &t, double Y[], double R[]) {
	double u = Y[2];
	double v = Y[3];

	double mod_v = sqrt(u * u + v * v);

	R[0] = u;
	R[1] = v;
	R[2] = -b * u * mod_v;
	R[3] = -1.0 - b * u * mod_v;
}

void shootingPlot(const double &thetaMin, const double &thetaMax,
                  const double &nTheta) {
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

void integrate(double y[], const double &theta, double xLast[], double yLast[],
               const int &order, std::ofstream &outFile) {
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

void integrate(double y[], const double &theta, double xLast[], double yLast[],
               const int &order) {
	std::ofstream dummyOutfile;
	integrate(y, theta, xLast, yLast, order, dummyOutfile);
}

void integrate(double y[], const double &theta, std::ofstream &outFile) {
	double dummyLast[2];
	integrate(y, theta, dummyLast, dummyLast, 1, outFile);
}

void integrate(double y[], const double &theta) {
	std::ofstream dummyOutfile;
	double dummyLast[2];
	integrate(y, theta, dummyLast, dummyLast, 1, dummyOutfile);
}

void integration(void (*RHSFunc)(const double &t, double *Y, double *RHS),
                 double y[], const int &nEq, double t, const double &dt,
                 const double &theta, double xLast[], double yLast[],
                 const int &order, const int &maxStep, std::ofstream &outFile) {
	outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
			<< "," << theta << endl;

	numIntegrations++;
	int stepCounter    = 0;
	bool exitCondition = false;
	while (stepCounter < maxStep && !exitCondition) {
		for (int i = 0; i < order - 1; i++) {
			xLast[i] = xLast[i + 1];
			yLast[i] = yLast[i + 1];
		}
		xLast[order - 1] = y[0];
		yLast[order - 1] = y[1];

		rk4Step(t, y, RHSFunc, dt, nEq);
		t += dt;
		stepCounter++;

		outFile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3]
				<< "," << theta << endl;

		if (xLast[0] < xTarg && y[0] > xTarg) exitCondition = true;
	}
}

double Residual(const double &theta) {
	double y[4];
	double xLast[64], yLast[64];
	if (gOrder > 64) throw exception("gOrder must be at most 64.");
	integrate(y, theta, xLast, yLast, gOrder);

	double xCurrent = y[0], yCurrent = y[1];

	return polInterp(xTarg, xLast, yLast, xCurrent, yCurrent, gOrder) - yTarg;
}
