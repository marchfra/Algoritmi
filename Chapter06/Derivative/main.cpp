#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/swap.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/quad.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/polynomials.hpp"
// #include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/root_finder.hpp"
#include "/Users/francescomarchisotti/Documents/Uni/Anno_4/Algoritmi/Libs/include/derivative.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

const double f(const double dX);
const double df(const double dX);

int main() {
	std::ofstream out;
	out.open("data.csv");
	if (!out) exit(1);
	out << "h,fd,bd,cd,hd\n";

	enum Method {None, ForwardDifference, BackwardDifference, CentralDifference, HigherDifference};
	std::string stringMethod[] = {"None", "ForwardDifference", "BackwardDifference", "CentralDifference", "HigherDifference"};

	double dX = 1.0;
	double dH = 0.5;
	const double dExact = df(dX);

	const int FW = 20, FS = 7, FP = 15;
	
	double dBest = std::numeric_limits<double>::max();
	double dBestDeriv, dBestH;
	Method BestMethod = None;

	while (dH > 1.0e-7) {
		double dDerivFD = ForwardDiff(f, dX, dH);
		double dErrorFD = dDerivFD - dExact;
		double dRelErrorFD = dErrorFD / dExact;
		if (fabs(dRelErrorFD) < dBest) {
			dBest = fabs(dRelErrorFD);
			dBestDeriv = dDerivFD;
			dBestH = dH;
			BestMethod = ForwardDifference;
		}

		double dDerivBD = BackwardDiff(f, dX, dH);
		double dErrorBD = dDerivBD - dExact;
		double dRelErrorBD = dErrorBD / dExact;
		if (fabs(dRelErrorBD) < dBest) {
			dBest = fabs(dRelErrorBD);
			dBestDeriv = dDerivBD;
			dBestH = dH;
			BestMethod = BackwardDifference;
		}

		double dDerivCD = CentralDiff(f, dX, dH);
		double dErrorCD = dDerivCD - dExact;
		double dRelErrorCD = dErrorCD / dExact;
		if (fabs(dRelErrorCD) < dBest) {
			dBest = fabs(dRelErrorCD);
			dBestDeriv = dDerivCD;
			dBestH = dH;
			BestMethod = CentralDifference;
		}

		double dDerivHD = HigherDiff(f, dX, dH);
		double dErrorHD = dDerivHD - dExact;
		double dRelErrorHD = dErrorHD / dExact;
		if (fabs(dRelErrorHD) < dBest) {
			dBest = fabs(dRelErrorHD);
			dBestDeriv = dDerivHD;
			dBestH = dH;
			BestMethod = HigherDifference;
		}

		cout.setf(std::ios::fixed);
		cout.precision(FP);
		cout << std::setw(FW) << dH
			 << std::setw(FS) << "FD"
			 << std::setw(FW) << dDerivFD
			 << std::setw(FW) << dRelErrorFD << "\n"
			 << std::setw(FW) << ""
			 << std::setw(FS) << "BD"
			 << std::setw(FW) << dDerivBD
			 << std::setw(FW) << dRelErrorBD << "\n"
			 << std::setw(FW) << ""
			 << std::setw(FS) << "CD"
			 << std::setw(FW) << dDerivCD
			 << std::setw(FW) << dRelErrorCD << "\n";

		dH /= 2.0;

		out << dH << "," << fabs(dErrorFD) << ","
						 << fabs(dErrorBD) << ","
						 << fabs(dErrorCD) << ","
						 << fabs(dErrorHD) << endl;
	}

	out.close();
	cout << "\n\n"
		 << "Best derivative estimate = " << dBestDeriv << endl
		 << "         corresponding h = " << dBestH << endl
		 << "             method used = " << stringMethod[BestMethod] << endl;

	return 0;
}

const double f(const double dX) {
	return sin(dX);
}

const double df(const double dX) {
	return cos(dX);
}
