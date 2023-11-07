#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

double gAlpha = 10.0;

const double f(const double t);

int main() {
	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);
	out << "t,f,df,ddf\n";

	int nN = 1000;
	double dT = 0.0;
	double dH = gAlpha / nN;

	while (dT <= gAlpha) {
		if (dT == 0.0)
			out << dT << "," << f(dT) << "," << ForwardDiff(f, dT, dH) << "," << ForwardDiff2(f, dT, dH) << endl;
		else 
			out << dT << "," << f(dT) << "," << CentralDiff(f, dT, dH) << "," << CentralDiff2(f, dT, dH) << endl;
		dT += dH;
	}

	out.close();
	return 0;
}

const double f(const double t) {
	if (t != 0.0)
		return gAlpha * t*t - t*t*t * (1 - exp(-gAlpha*gAlpha / t));
	else if (t < 1.0e-3)
		return t*t * (gAlpha - t);
	else
		return 0;
}
