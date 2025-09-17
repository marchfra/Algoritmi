#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "../../Libs/include/swap.hpp"
// #include "../../Libs/include/quad.hpp"
// #include "../../Libs/include/polynomials.hpp"
// #include "../../Libs/include/root_finder.hpp"
// #include "../../Libs/include/derivative.hpp"
// #include "../../Libs/include/ode_solver.hpp"
// #include "../../Libs/include/lin_alg.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

int main() {
	std::ofstream out;
	out.open("../data/data.csv");
	if (!out) exit(5);

	cout << "Hello World!" << endl;

	out.close();
	return 0;
}
