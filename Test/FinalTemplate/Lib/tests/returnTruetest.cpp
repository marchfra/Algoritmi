#include "test_config.hpp"
#include "../include/mycode.hpp"

TEST_CASE("returnTrue test") {
	CHECK(returnTrue());
	CHECK_FALSE(returnTrue(5));
}
