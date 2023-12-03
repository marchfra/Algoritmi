#include "../include/quad.hpp"

#include "../include/debug.hpp"

double rectangularQuad(double (*F)(const double& x), const double& a,
                       const double& b, const int& n) {
	double sum     = 0.0;
	const double h = (b - a) / n;  // Interval width
	double xi      = a;

	for (int i = 0; i < n; i++) {  // Rectangular rule
		sum += F(xi);
		xi += h;
	}

	return sum * h;  // Finalise the calculation and return
}

double midpointQuad(double (*F)(const double& x), const double& a,
                    const double& b, const int& n) {
	return gaussLegendreQuad(
		F, a, b, n,
		1);  // Gauss-Legendre rule converges to midpoint rule for Ng = 1
}

double trapezoidalQuad(double (*F)(const double& x), const double& a,
                       const double& b, const int& n) {
	double sum =
		0.5 *
		(F(a) + F(b));  // Take care of the first and last term of the summation
	const double h = (b - a) / n;  // Interval width
	double xi      = a;

	for (int i = 1; i < n; i++) {  // Trapezoidal rule
		xi += h;
		sum += F(xi);
	}

	return sum * h;  // Finalise the calculation and return
}

double simpsonQuad(double (*F)(const double& x), const double& a,
                   const double& b, const int& n) {
	// n must be even!
	if (n % 2 != 0) {
		throw std::invalid_argument(
			"simpsonQuad(): Invalid argument: n must be even.");
	}

	const double h = (b - a) / n;  // Interval width
	double sum =
		F(a) + F(b);  // Take care of the first and last term of the summation
	double xi = a;
	double w  = 4.0;  // Define the second weight to be 4

	for (int i = 1; i < n; i++) {
		xi += h;
		sum += F(xi) * w;
		w = 6.0 - w;  // Alternate between w = 4 and w = 2
	}

	return sum * h / 3;  // Finalise the calculation and return
}

void setLegendreWeightsAndRoots(double weights[], double roots[],
                                const int& Ng) {
	switch (Ng) {
	case 1:
		roots[0]   = 0;
		weights[0] = 2;
		break;
	case 2:
		roots[0]   = sqrt(1.0 / 3.0);
		weights[0] = 1;
		roots[1]   = -sqrt(1.0 / 3.0);
		weights[1] = 1;
		break;
	case 3:
		roots[0]   = 0;
		weights[0] = 8.0 / 9.0;
		roots[1]   = sqrt(3.0 / 5.0);
		weights[1] = 5.0 / 9.0;
		roots[2]   = -sqrt(3.0 / 5.0);
		weights[2] = 5.0 / 9.0;
		break;
	case 4:
		roots[0]   = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		weights[0] = (18 + sqrt(30)) / 36;
		roots[1]   = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		weights[1] = (18 + sqrt(30)) / 36;
		roots[2]   = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		weights[2] = (18 - sqrt(30)) / 36;
		roots[3]   = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		weights[3] = (18 - sqrt(30)) / 36;
		break;
	case 5:
		roots[0]   = 0;
		weights[0] = 128.0 / 225.0;
		roots[1]   = 1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));
		weights[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		roots[2]   = -1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));
		weights[2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		roots[3]   = 1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));
		weights[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		roots[4]   = -1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));
		weights[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		break;
	default:
		throw std::invalid_argument("setLegendreWeightsAndRoots(): Invalid "
		                            "argument: Ng must be in [1, 5].");
	}
}

double gaussLegendreQuad(double (*F)(const double& x), const double& a,
                         const double& b, const int N, const int Ng) {
	double weight[8], root[8];

	setLegendreWeightsAndRoots(weight, root, Ng);

	const double h = (b - a) / N;  // Interval width
	double sum     = 0.0;

	double xc, sumj;

	for (int i = 0; i < N; i++) {  // Loop over intervals
		xc   = a + (i + 0.5) * h;  // Define the centre of the interval
		sumj = 0.0;                // Initialize sum for this interval
		for (int j = 0; j < Ng;
		     j++) {  // Loop over gaussian points in the interval
			sumj += weight[j] *
			        F(0.5 * h * root[j] + xc);  // and apply gaussian rule
		}
		sum += sumj;  // Add partial sum to total integral
	}

	return 0.5 * h * sum;  // Finalize calculation and return
}

double gaussLegendreQuad2D(double (*F)(const double& x, const double& y),
                           const double& xa, const double& xb, const double& ya,
                           const double& yb, const int N, const int Ng) {
	double weight[8], root[8];

	setLegendreWeightsAndRoots(weight, root, Ng);

	const double xh = (xb - xa) / N;  // Interval width along x
	const double yh = (yb - ya) / N;  // Interval width along y
	double sum      = 0.0;

	double xc, yc, xsum, ysum;

	for (int i = 0; i < N; i++) {      // Loop over x intervals
		for (int j = 0; j < N; j++) {  // Loop over y intervals
			xc = xa + (i + 0.5) * xh;  // Define the centre of the
			yc = ya + (j + 0.5) * yh;  // x and y intervals

			ysum = 0.0;  // Initialize sum over x for this interval
			for (int jk = 0; jk < Ng; jk++) {  // Loop over y gaussian points
				xsum = 0.0;  // Initialize sum over y for this interval
				for (int ik = 0; ik < Ng;
				     ik++) {  // Loop over x gaussian points and apply gaussian
					          // rule
					xsum += weight[ik] * F(0.5 * xh * root[ik] + xc,
					                       0.5 * yh * root[jk] + yc);
				}
				ysum += weight[jk] * xsum;  // Add weighted xsum to y sum
			}
			sum += ysum;  // Add partial sum to total integral
		}
	}

	return 0.25 * xh * yh * sum;  // Finalise calculation and return
}
