#include <cmath>

#include "test_config.hpp"
#include "../include/derivative.hpp"

double func1(const double& x);
double dfunc1(const double& x);
double ddfunc1(const double& x);

double func2(const double& x);
double dfunc2(const double& x);
double ddfunc2(const double& x);

TEST_CASE("testing first derivative functions") {
	double x[] = {0.35264398388043505, 0.3223567365282518, 0.5341485343520432, 0.0699975902021307, 0.7386715788205122};
	int n = static_cast<int>(sizeof(x) / sizeof(x[0]));
	SUBCASE("testing on func1") {
		double h = 0.1;

		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])).epsilon(0.1));
			CHECK(backwardDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])).epsilon(0.1));
			CHECK(centralDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])).epsilon(0.05));
			CHECK(higherDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])));
		}

		h = 1.0e-5;
		
		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])));
			CHECK(backwardDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])));
			CHECK(centralDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])));
			CHECK(higherDiff(func1, x[i], h) == doctest::Approx(dfunc1(x[i])));
		}
	}

	SUBCASE("testing on func2") {
		double h = 0.1;

		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])).epsilon(0.1));
			CHECK(backwardDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])).epsilon(0.1));
			CHECK(centralDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])).epsilon(0.05));
			CHECK(higherDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])));
		}

		h = 1.0e-5;
		
		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])));
			CHECK(backwardDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])));
			CHECK(centralDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])));
			CHECK(higherDiff(func2, x[i], h) == doctest::Approx(dfunc2(x[i])));
		}
	}
}

TEST_CASE("testing second derivative functions") {
	double x[] = {0.35264398388043505, 0.3223567365282518, 0.5341485343520432, 0.0699975902021307, 0.7386715788205122};
	int n = static_cast<int>(sizeof(x) / sizeof(x[0]));
	SUBCASE("testing on func1") {
		double h = 0.1;

		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])).epsilon(0.1));
			CHECK(backwardDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])).epsilon(0.1));
			CHECK(centralDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])).epsilon(0.05));
		}

		h = 1.0e-5;
		
		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])));
			CHECK(backwardDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])));
			CHECK(centralDiff2(func1, x[i], h) == doctest::Approx(ddfunc1(x[i])));
		}
	}

	SUBCASE("testing on func2") {
		double h = 0.1;

		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])).epsilon(0.1));
			CHECK(backwardDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])).epsilon(0.1));
			CHECK(centralDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])).epsilon(0.05));
		}

		h = 1.0e-5;
		
		for (int i = 0; i < n; i++) {
			CHECK(forwardDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])));
			CHECK(backwardDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])));
			CHECK(centralDiff2(func2, x[i], h) == doctest::Approx(ddfunc2(x[i])));
		}
	}
}

double func1(const double& x) {
	return sin(x);
}

double dfunc1(const double& x) {
	return cos(x);
}

double ddfunc1(const double& x) {
	return -sin(x);
}

double func2(const double& x) {
	return x * exp(-x);
}

double dfunc2(const double& x) {
	return exp(-x) * (1.0 - x);
}

double ddfunc2(const double& x) {
	return exp(-x) * (x - 2.0);
}
