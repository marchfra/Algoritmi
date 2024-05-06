#include <cmath>
#include <exception>

#include "test_config.hpp"
#include "../include/root_finder_param.hpp"
#include "../include/polynomials.hpp"

double func1(const double& x, const double& k);
double dfunc1(const double& x, const double& k);

double func2(const double& x, const double& k);
double dfunc2(const double& x, const double& k);

double func3(const double& x, const double& k);
double dfunc3(const double& x, const double& k);

double func4(const double& x, const double& k);
double dfunc4(const double& x, const double& k);

double func5(const double& x, const double& k);
double dfunc5_small(const double& x, const double& k);
double dfunc5(const double& x, const double& k);

TEST_CASE("testing bisection function") {
	double xtol = 1.0e-7;
	double xa = -1.0, xb = 1.0;

	double root = 0.0;

	SUBCASE("testing exceptions") {
		CHECK_THROWS_WITH_AS(bisection(func1, 1.0, xa, -xb, xtol, root),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;

		bisection(func1, 1.0, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		bisection(func1, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 26);

		double ftol = 0.5;
		bisection(func1, 1.0, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 26);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		bisection(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 30);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		bisection(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		bisection(func3, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 26);
	}
}

TEST_CASE("testing falsePosition function") {
	double xtol = 1.0e-7;
	double xa = -1.0, xb = 1.0;

	double root = 0.0;

	SUBCASE("testing exceptions") {
		CHECK_THROWS_WITH_AS(falsePosition(func1, 1.0, xa, -xb, xtol, root),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;

		falsePosition(func1, 1.0, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		falsePosition(func1, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 15);

		double ftol = 0.5;
		falsePosition(func1, 1.0, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 15);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		falsePosition(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 80);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		falsePosition(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 22);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		falsePosition(func3, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 55);
	}
}

TEST_CASE("testing secant function") {
	double xtol;
	double xa, xb;

	double root = 0.0;


	SUBCASE("testing exceptions") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;

		CHECK_THROWS_WITH_AS(secant(func3, 1.0, xa, xb, xtol, root),
							 "Maximum number of steps exceeded.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671438e-01;

		secant(func1, 1.0, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		secant(func1, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 6);

		double ftol = 0.5;
		secant(func1, 1.0, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 6);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		secant(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 12);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		secant(func2, 1.0, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 10);
	}
}

TEST_CASE("testing newton function") {
	double xtol;
	double xa, xb;
	double dftol = 1.0e-3;

	double root = 0.0;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;

		newton(func1, dfunc1, 1.0, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		newton(func1, dfunc1, 1.0, xa, xb, xtol, dftol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 5);

		double ftol = 0.5;
		newton(func1, dfunc1, 1.0, xa, xb, xtol, ftol, dftol, root, nTry);
		CHECK(nTry < 5);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		newton(func2, dfunc2, 1.0, xa, xb, xtol, dftol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 6);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		newton(func2, dfunc2, 1.0, xa, xb, xtol, dftol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		newton(func3, dfunc3, 1.0, xa, xb, xtol, dftol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 8);
	}

	SUBCASE("testing derivative tolerance") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		double ftol = -1.0;
		const double expected = 1.0;

		CHECK_THROWS_WITH_AS(newton(func5, dfunc5_small, 1.0, xa, xb, xtol, root),
							 "Derivative too small.",
							 std::runtime_error);

		int nTry = -1;
		newton(func5, dfunc5, 1.0, xa, xb, xtol, ftol, dftol, root, nTry);
		CHECK(root == expected);
	}
}

TEST_CASE("testing bracket function") {
	const int N = 10;
	const double xa = -10.0, xb = 10.0;

	double xL[8], xR[8];
	int nRoots;

	const double xLExpected[] = {-10.0, -8.0, -4.0, 0.0, 2.0};
	const double xRExpected[] = {-8.0, -6.0, -2.0, 2.0, 4.0};
	const int nRootsExpected = 5;

	bracket(func4, 1.0, xa, xb, xL, xR, N, nRoots);

	REQUIRE(nRoots == nRootsExpected);
	for (int i = 0; i < nRoots; i++) {
		CHECK(xL[i] == doctest::Approx(xLExpected[i]));
		CHECK(xR[i] == doctest::Approx(xRExpected[i]));
	}
}

TEST_CASE("testing findRoots function") {
	int N = 10;
	double xa = -10.0, xb = 10.0;
	const double tol = 1.0e-7;
	double roots[8];
	int nRoots;

	const int nRootsExpected = 5;
	const double rootsExpected[] = {-8.716925e+00, -6.889594e+00, -2.968485e+00, 4.361680e-01, 2.183971e+00};

	SUBCASE("testing exceptions") {
		N = 1;
		xa = 0.0; xb = 2.0;
		CHECK_THROWS_WITH_AS(findRoots(func3, 1.0, xa, xb, tol, roots, nRoots, N, "secant"),
							 "Maximum number of steps exceeded.",
							 std::runtime_error);

		xa = 5.0; xb = 10.0;
		CHECK_THROWS_WITH_AS(findRoots(func1, 1.0, xa, xb, tol, roots, nRoots, N),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);

		xa = -10.0; xb = 10.0;
		CHECK_THROWS_WITH_AS(findRoots(func1, 1.0, xa, xb, tol, roots, nRoots, N, "gobble"),
							 "Invalid method argument.",
							 std::invalid_argument);
	}

	SUBCASE("testing newton implementation") {
		findRoots(func4, dfunc4, 1.0, xa, xb, tol, roots, nRoots, N);

		CHECK(nRoots == nRootsExpected);
		for (int i = 0; i < nRoots; i++) {
			CHECK(roots[i] == doctest::Approx(rootsExpected[i]));
		}
	}

	SUBCASE("testing overloaded") {
		findRoots(func4, 1.0, xa, xb, tol, roots, nRoots, N);

		CHECK(nRoots == nRootsExpected);
		for (int i = 0; i < nRoots; i++) {
			CHECK(roots[i] == doctest::Approx(rootsExpected[i]));
		}
	}
}

double func1(const double& x, const double& k) {
	return exp(-x) - x;
}

double dfunc1(const double& x, const double& k) {
	return -exp(-x) - 1.0;
}

double func2(const double& x, const double& k) {
	const double a[] = {5.0, 1.0, -3.0, 1.0};
	const int degree = 3;
	return hornerPol(x, a, degree);
}

double dfunc2(const double& x, const double& k) {
	const double a[] = {5.0, 1.0, -3.0, 1.0};
	const int degree = 3;
	double dpdx;
	hornerPol(x, a, degree, dpdx);
	return dpdx;
}

double func3(const double& x, const double& k) {
	return exp(1 / (x + 0.5)) - (3 + 2 * x) / (1 + x);
}

double dfunc3(const double& x, const double& k) {
	return 1 / ((1 + x)*(1 + x)) - exp(1 / (x + 0.5)) / ((x + 0.5)*(x + 0.5));
}

double func4(const double& x, const double& k) {
	double bracket = (0.1 * x)*(0.1 * x) + 0.2 * x + 1.0 / 3;
	return sin(x) - bracket;
}

double dfunc4(const double& x, const double& k) {
	double bracket = 2.0 * (0.1 * x) * 0.1 + 0.2;
	return cos(x) - bracket;
}

double func5(const double& x, const double& k) {
	return x*x - 1.0;
}

double dfunc5_small(const double& x, const double& k) {
	return 1.0e-8;
}

double dfunc5(const double& x, const double& k) {
	return 2.0 * x;
}
