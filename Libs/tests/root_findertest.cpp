#include <cmath>

#include "test_config.hpp"
#include "../include/root_finder.hpp"
#include "../include/polynomials.hpp"

double func1(const double& x);
double dfunc1(const double& x);

double func2(const double& x);
double dfunc2(const double& x);

double func3(const double& x);
double dfunc3(const double& x);

double func4(const double& x);
double dfunc4(const double& x);

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

TEST_CASE("testing bracke function") {
	const int N = 10;
	const double xa = -10.0, xb = 10.0;

	double xL[8], xR[8];
	int nRoots;

	const double xLExpected[] = {-10.0, -8.0, -4.0, 0.0, 2.0};
	const double xRExpected[] = {-8.0, -6.0, -2.0, 2.0, 4.0};
	const int nRootsExpected = 5;

	bracket(func4, xa, xb, xL, xR, N, nRoots);

	REQUIRE(nRoots == nRootsExpected);
	for (int i = 0; i < nRoots; i++) {
		CHECK(xL[i] == doctest::Approx(xLExpected[i]));
		CHECK(xR[i] == doctest::Approx(xRExpected[i]));
	}
}

TEST_CASE("testing findRoots function") {
	const int N = 10;
	const double xa = -10.0, xb = 10.0;
	const double tol = 1.0e-7;
	double roots[8];
	int nRoots;

	const int nRootsExpected = 5;
	const double rootsExpected[] = {-8.716925e+00, -6.889594e+00, -2.968485e+00, 4.361680e-01, 2.183971e+00};

	SUBCASE("testing newton implementation") {
		findRoots(func4, dfunc4, xa, xb, tol, roots, nRoots, N);

		REQUIRE(nRoots == nRootsExpected);
		for (int i = 0; i < nRoots; i++) {
			CHECK(roots[i] == doctest::Approx(rootsExpected[i]));
		}
	}

	SUBCASE("testing overloaded") {
		findRoots(func4, xa, xb, tol, roots, nRoots, N);

		REQUIRE(nRoots == nRootsExpected);
		for (int i = 0; i < nRoots; i++) {
			CHECK(roots[i] == doctest::Approx(rootsExpected[i]));
		}
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

double func4(const double& x) {
	double bracket = (0.1 * x)*(0.1 * x) + 0.2 * x + 1.0 / 3;
	return sin(x) - bracket;
}

double dfunc4(const double& x) {
	double bracket = 2.0 * (0.1 * x) * 0.1 + 0.2;
	return cos(x) - bracket;
}
