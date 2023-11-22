// =============================================================================
// =============================================================================
// ================== findRoots() AND bracket() AREN'T TESTED ==================
// =============================================================================
// =============================================================================

#include <cmath>
// #include <limits>

#include "test_config.hpp"
#include "../include/root_finder.hpp"
#include "../include/polynomials.hpp"

double func1(const double& x);
double dfunc1(const double& x);

double func2(const double& x);
double dfunc2(const double& x);

double func3(const double& x);
double dfunc3(const double& x);

TEST_CASE("testing bisection function") {
	double xtol;
	double xa, xb;

	double root = 0.0;
	int flag = -1;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		flag = bisection(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);

		int nTry = 0;
		flag = bisection(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 26);

		double ftol = 0.5;
		flag = bisection(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(flag == 0);
		CHECK(nTry < 26);

		flag = bisection(func1, xa, -xb, xtol, root);
		CHECK(flag == 2);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		flag = bisection(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 30);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		flag = bisection(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		flag = bisection(func3, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 26);
	}
}

TEST_CASE("testing falsePosition function") {
	double xtol;
	double xa, xb;

	double root = 0.0;
	int flag = -1;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		flag = falsePosition(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);

		int nTry = 0;
		flag = falsePosition(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 15);

		double ftol = 0.5;
		flag = falsePosition(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(flag == 0);
		CHECK(nTry < 15);

		flag = falsePosition(func1, xa, -xb, xtol, root);
		CHECK(flag == 2);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		flag = falsePosition(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 80);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		flag = falsePosition(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 22);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		flag = falsePosition(func3, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 55);
	}
}

TEST_CASE("testing secant function") {
	double xtol;
	double xa, xb;

	double root = 0.0;
	int flag = -1;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671438e-01;
		
		flag = secant(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);

		int nTry = 0;
		flag = secant(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 6);

		double ftol = 0.5;
		flag = secant(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(flag == 0);
		CHECK(nTry < 6);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		flag = secant(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 12);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		flag = secant(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 10);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;

		int nTry = 0;
		flag = secant(func3, xa, xb, xtol, root, nTry);
		CHECK(flag == 1);
	}
}

TEST_CASE("testing newton function") {
	double xtol;
	double xa, xb;

	double root = 0.0;
	int flag = -1;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		flag = newton(func1, dfunc1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);

		int nTry = 0;
		flag = newton(func1, dfunc1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 5);

		double ftol = 0.5;
		flag = newton(func1, dfunc1, xa, xb, xtol, ftol, root, nTry);
		CHECK(flag == 0);
		CHECK(nTry < 5);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		flag = newton(func2, dfunc2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 6);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		flag = newton(func2, dfunc2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		flag = newton(func3, dfunc3, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(flag == 0);
		CHECK(nTry == 8);
	}
}

double func1(const double& x) {
	return exp(-x) - x;
}

double dfunc1(const double& x) {
	return -exp(-x) - 1.0;
}

double func2(const double& x) {
	const double a[] = {5.0, 1.0, -3.0, 1.0};
	const int degree = 3;
	return hornerPol(x, a, degree);
}

double dfunc2(const double& x) {
	const double a[] = {5.0, 1.0, -3.0, 1.0};
	const int degree = 3;
	double dpdx;
	hornerPol(x, a, degree, dpdx);
	return dpdx;
}

double func3(const double& x) {
	return exp(1 / (x + 0.5)) - (3 + 2 * x) / (1 + x);
}

double dfunc3(const double& x) {
	return 1 / ((1 + x)*(1 + x)) - exp(1 / (x + 0.5)) / ((x + 0.5)*(x + 0.5));
}
