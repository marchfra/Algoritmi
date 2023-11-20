#include <iostream>

#include "mycode.hpp"

int main(int argc, char **argv) {
	if (argc > 0) {
		std::cout << argv[1] << std::endl;
		std::cout << reverse(argv[1]) << std::endl;
	}

	return 0;
}
