#include "test_config.hpp"
#include "../include/swap.hpp"

TEST_CASE("testing swap function") {
	float a = 3.14, b = 1.62;
	swap(a, b);
	CHECK(a == doctest::Approx(1.62));
	CHECK(b == doctest::Approx(3.14));
}

TEST_CASE("testing swap function 2") {
	double a = 3.14, b = 1.62;
	swap(a, b);
	CHECK(a == doctest::Approx(1.62));
	CHECK(b == doctest::Approx(3.14));
}
