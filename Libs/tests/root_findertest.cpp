#include <cmath>
#include <exception>

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
	double xtol = 1.0e-7;
	double xa = -1.0, xb = 1.0;

	double root = 0.0;

	SUBCASE("testing exceptions") {
		CHECK_THROWS_WITH_AS(bisection(func1, xa, -xb, xtol, root),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		bisection(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		bisection(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 26);

		double ftol = 0.5;
		bisection(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 26);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		bisection(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 30);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		bisection(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		bisection(func3, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 26);
	}
}

TEST_CASE("testing falsePosition function") {
	double xtol = 1.0e-7;
	double xa = -1.0, xb = 1.0;

	double root = 0.0;

	SUBCASE("testing exceptions") {
		CHECK_THROWS_WITH_AS(falsePosition(func1, xa, -xb, xtol, root),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		falsePosition(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		falsePosition(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 15);

		double ftol = 0.5;
		falsePosition(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 15);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		falsePosition(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 80);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		falsePosition(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 22);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		falsePosition(func3, xa, xb, xtol, root, nTry);
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

		CHECK_THROWS_WITH_AS(secant(func3, xa, xb, xtol, root),
							 "Maximum number of steps exceeded.",
							 std::runtime_error);
	}

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671438e-01;
		
		secant(func1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		secant(func1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 6);

		double ftol = 0.5;
		secant(func1, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 6);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		secant(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 12);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		secant(func2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 10);
	}
}

TEST_CASE("testing newton function") {
	double xtol;
	double xa, xb;

	double root = 0.0;

	SUBCASE("roots of func1") {
		xtol = 1.0e-7;
		xa = -1.0; xb = 1.0;
		const double expected = 5.671433e-01;
		
		newton(func1, dfunc1, xa, xb, xtol, root);
		CHECK(root == doctest::Approx(expected));

		int nTry = 0;
		newton(func1, dfunc1, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 5);

		double ftol = 0.5;
		newton(func1, dfunc1, xa, xb, xtol, ftol, root, nTry);
		CHECK(nTry < 5);
	}

	SUBCASE("roots of func2") {
		xtol = 1.0e-8;
		xa = -5.0; xb = 0.0;
		const double expected = -1.0;

		int nTry = 0;
		newton(func2, dfunc2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 6);

		xa = -2.0; xb = 0.0;
		nTry = 0;
		newton(func2, dfunc2, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 1);
	}

	SUBCASE("roots of func3") {
		xtol = 1.0e-7;
		xa = 0.0; xb = 2.0;
		const double expected = 5.235934e-01;

		int nTry = 0;
		newton(func3, dfunc3, xa, xb, xtol, root, nTry);
		CHECK(root == doctest::Approx(expected));
		CHECK(nTry == 8);
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

	bracket(func4, xa, xb, xL, xR, N, nRoots);

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
		CHECK_THROWS_WITH_AS(findRoots(func3, xa, xb, tol, roots, nRoots, N, "secant"),
							 "Maximum number of steps exceeded.",
							 std::runtime_error);

		xa = 5.0; xb = 10.0;
		CHECK_THROWS_WITH_AS(findRoots(func1, xa, xb, tol, roots, nRoots, N),
							 "The supplied interval does not contain any roots.",
							 std::runtime_error);

		xa = -10.0; xb = 10.0;
		CHECK_THROWS_WITH_AS(findRoots(func1, xa, xb, tol, roots, nRoots, N, "gobble"),
							 "Invalid method argument.",
							 std::invalid_argument);
	}

	SUBCASE("testing newton implementation") {
		findRoots(func4, dfunc4, xa, xb, tol, roots, nRoots, N);

		CHECK(nRoots == nRootsExpected);
		for (int i = 0; i < nRoots; i++) {
			CHECK(roots[i] == doctest::Approx(rootsExpected[i]));
		}
	}

	SUBCASE("testing overloaded") {
		findRoots(func4, xa, xb, tol, roots, nRoots, N);

		CHECK(nRoots == nRootsExpected);
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
