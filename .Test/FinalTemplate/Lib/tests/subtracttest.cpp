#include "test_config.hpp"
#include "../include/subtract.hpp"

TEST_CASE("testing subtract function") {
	int a = 13, b = 4;
	CHECK(subtract(a, b) == 9);
}
