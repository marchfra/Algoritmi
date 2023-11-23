#include "test_config.hpp"
#include "../include/polynomials.hpp"

TEST_CASE("testing hornerPol function") {
	SUBCASE("testing with f(x) = x^3 - 3x^2 + x + 5") {
		const double a[] = {5.0, 1.0, -3.0, 1.0};
		const int degree = 3;
		double dpdx;

		double x = -1.0;
		double expectedPol = 0.0;
		double expectedDPol = 10.0;
		double pol = hornerPol(x, a, degree, dpdx);
		CHECK(pol == doctest::Approx(expectedPol));
		CHECK(dpdx == doctest::Approx(expectedDPol));

		x = 0.0;
		expectedPol = 5.0;
		expectedDPol = 1.0;
		pol = hornerPol(x, a, degree, dpdx);
		CHECK(pol == doctest::Approx(expectedPol));
		CHECK(dpdx == doctest::Approx(expectedDPol));

		x = 5.0;
		expectedPol = 60.0;
		expectedDPol = 46.0;
		pol = hornerPol(x, a, degree, dpdx);
		CHECK(pol == doctest::Approx(expectedPol));
		CHECK(dpdx == doctest::Approx(expectedDPol));

		x = -2.7;
		expectedPol = -39.253;
		expectedDPol = 39.07;
		pol = hornerPol(x, a, degree, dpdx);
		CHECK(pol == doctest::Approx(expectedPol));
		CHECK(dpdx == doctest::Approx(expectedDPol));

		x = 3.14;
		expectedPol = 9.520344;
		pol = hornerPol(x, a, degree);
		CHECK(pol == doctest::Approx(expectedPol));
	}
}
