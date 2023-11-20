#include "test_config.hpp"
#include "../include/add.hpp"

TEST_CASE("testing add function") {
	int a = 13, b = 4;
	CHECK(add(a, b) == 17);
}
