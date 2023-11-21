#include <cmath>
#include <limits>

#include "test_config.hpp"
#include "../include/quad.hpp"

double func1(const double& x);
double func2(const double& x);
double func3(const double& x, const double& y);

TEST_CASE("testing rectangularQuad function") {
	int N = 4;
	const double xa = 0.0, xb = 1.0;

	SUBCASE("fixed subinterval number") {
		const double expected = 7.144244988813e-01;
		
		CHECK(rectangularQuad(func1, xa, xb, N) == doctest::Approx(expected));
	}

	SUBCASE("convergence") {
		const double expected = 6.321302042719e-01;
		const int expectedN = 32'768;
		const double tol = 1.0e-5;
		double err = std::numeric_limits<double>::max();;
		double prevI = std::numeric_limits<double>::max();;

		while (err > tol) {
			const double I = rectangularQuad(func1, xa, xb, N);
			err = fabs(I - prevI);
			prevI = I;
			N *= 2;
		}
		N /= 2;

		CHECK(prevI == doctest::Approx(expected));
		CHECK(N == expectedN);
	}
}

TEST_CASE("testing trapezoidalQuad function") {
	int N = 4;
	const double xa = 0.0, xb = 1.0;

	SUBCASE("fixed subinterval number") {
		const double expected = 6.354094290277e-01;
		
		CHECK(trapezoidalQuad(func1, xa, xb, N) == doctest::Approx(expected));
	}

	SUBCASE("convergence") {
		const double expected = 6.321237739567e-01;
		const int expectedN = 128;
		const double tol = 1.0e-5;
		double err = std::numeric_limits<double>::max();;
		double prevI = std::numeric_limits<double>::max();;

		while (err > tol) {
			const double I = trapezoidalQuad(func1, xa, xb, N);
			err = fabs(I - prevI);
			prevI = I;
			N *= 2;
		}
		N /= 2;

		CHECK(prevI == doctest::Approx(expected));
		CHECK(N == expectedN);
	}
}

TEST_CASE("testing simpsonQuad function") {
	int N = 4;
	double xa = 0.0, xb = 1.0;

	SUBCASE("parity of the number of subintervals") {
		N = 3;

		CHECK_THROWS_WITH_AS(simpsonQuad(func1, xa, xb, N),
							 "simpsonQuad(): Invalid argument: n must be even.",
							 std::invalid_argument);

		N = 2;
		CHECK_NOTHROW(simpsonQuad(func1, xa, xb, N));
	}

	SUBCASE("fixed subinterval number") {
		const double expected = 6.321341753205e-01;
		
		CHECK(simpsonQuad(func1, xa, xb, N) == doctest::Approx(expected));
	}

	SUBCASE("convergence") {
		const double expected = 6.321206123892e-01;
		const int expectedN = 16;
		const double tol = 1.0e-5;
		double err = std::numeric_limits<double>::max();;
		double prevI = std::numeric_limits<double>::max();;

		while (err > tol) {
			const double I = simpsonQuad(func1, xa, xb, N);
			err = fabs(I - prevI);
			prevI = I;
			N *= 2;
		}
		N /= 2;

		CHECK(prevI == doctest::Approx(expected));
		CHECK(N == expectedN);
	}

	SUBCASE("test on function 2: sqrt(1 + x)") {
		xa = 0.0; xb = 3.0;
		N = 2;
		const double expected = 4.662277660168e+00;

		CHECK(simpsonQuad(func2, xa, xb, N) == doctest::Approx(expected));
	}
}

TEST_CASE("testing setLegendreWeightsAndRoots function") {
	double weights[8], roots[8];

	SUBCASE("number of gaussian points greater than max") {
		CHECK_THROWS_WITH_AS(setLegendreWeightsAndRoots(weights, roots, 8), 
							 "setLegendreWeightsAndRoots(): Invalid argument: Ng must be in [1, 5].",
							 std::invalid_argument);
	}

	SUBCASE("normal operation") {
		const double expectedRoot1 = 0.0;
		const double expectedRoot2 = 0.577350269189626;
		const double expectedRoot3[] = {0.0, 0.774596669241483};

		const double expectedWeight1 = 2.0;
		const double expectedWeight2 = 1.0;
		const double expectedWeight3[] = {0.888888888888889, 0.555555555555556};

		for (int Ng = 1; Ng <= 3; Ng++) {
			setLegendreWeightsAndRoots(weights, roots, Ng);
			switch (Ng) {
			case 1:
				CHECK(roots[0] == doctest::Approx(expectedRoot1));
				CHECK(weights[0] == doctest::Approx(expectedWeight1));
				break;
			case 2:
				CHECK(roots[0] == doctest::Approx(expectedRoot2));
				CHECK(weights[0] == doctest::Approx(expectedWeight2));
				break;
			case 3:
				for (int i = 0; i < 2; i++) {
					CHECK(roots[i] == doctest::Approx(expectedRoot3[i]));
					CHECK(weights[i] == doctest::Approx(expectedWeight3[i]));
				}
				break;
			default:
				exit(1);
			}
		}
	}
}

TEST_CASE("testing gaussLegendreQuad function") {
	const int N = 1;
	const int Ng = 3;
	const double xa = 0.0, xb = 3.0;
	
	const double expected = 4.666829051581e+00;

	SUBCASE("test on function 2: sqrt(1 + x)") {
		CHECK(gaussLegendreQuad(func2, xa, xb, N, Ng) == doctest::Approx(expected));
	}

	SUBCASE("test default gaussian points number") {
		CHECK(gaussLegendreQuad(func2, xa, xb, N) == doctest::Approx(expected));
	}

	SUBCASE("test default subintervals number") {
		CHECK(gaussLegendreQuad(func2, xa, xb) == doctest::Approx(expected));
	}
}

TEST_CASE("testing gaussLegendreQuad2D function") {
	const double xa = -1.0, xb = 1.0;
	const double ya = -1.0, yb = 1.0;
	const int N = 1;
	const int Ng = 3;

	const double expected = 412.0 / 45.0;

	SUBCASE("test on function 2: sqrt(1 + x)") {
		CHECK(gaussLegendreQuad2D(func3, xa, xb, ya, yb, N, Ng) == doctest::Approx(expected));
	}

	SUBCASE("test default gaussian points number") {
		CHECK(gaussLegendreQuad2D(func3, xa, xb, ya, yb, N) == doctest::Approx(expected));
	}

	SUBCASE("test default subintervals number") {
		CHECK(gaussLegendreQuad2D(func3, xa, xb, ya, yb) == doctest::Approx(expected));
	}
}

double func1(const double& x) {
	return exp(-x);
}

double func2(const double& x) {
	return sqrt(1 + x);
}

double func3(const double& x, const double& y) {
	return x*x*x*x * y*y + 2.0 * x*x * y*y - x*x * y + 2.0;
}
