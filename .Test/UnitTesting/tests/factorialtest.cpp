#include "test_config.hpp"
#include "../include/mycode.hpp"

TEST_CASE("testing the factorial function") {
    CHECK(factorial(0) == 1);
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
    CHECK_THROWS(factorial(-2));
}
