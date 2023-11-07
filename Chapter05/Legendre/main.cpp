#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
#include "GaussLegendre.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

int g_nEval = 0;

double legendre(const double& x, const int& n);

int main() {
	cout << "Hello World!" << endl;
	cout.precision(5);

	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);

	out << "n,nEval" << endl;

	for (int n = 0; n < 30; n++) {	
		cout << std::setw(9) << legendre(0.1, n) << "   n = " << std::setw(2) << n << "   nEval = " << g_nEval << endl;
		out << n << "," << g_nEval << endl;
		g_nEval = 0;
	}

	out.close();

	return 0;
}

double legendre(const double& x, const int& n) {
	g_nEval++;
	switch (n) {
	case 1:
		return x;
	case 0:
		return 1;
	default:
		return static_cast<double>(2 * n - 1) / n * x * legendre(x, n - 1) - static_cast<double>(n - 1) / n * legendre(x, n - 2);
	}
}
