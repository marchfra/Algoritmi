#include "../include/debug.hpp"
#include "../include/derivative.hpp"

double ForwardDiff(const double (*f)(const double), const double dX, const double dH) {
	double fxph = f(dX + dH);
	double fx   = f(dX);
	
	return (fxph - fx) / dH;
}

double BackwardDiff(const double (*f)(const double), const double dX, const double dH) {
	double fx   = f(dX);
	double fxmh = f(dX - dH);
	
	return (fx - fxmh) / dH;
}

double CentralDiff(const double (*f)(const double), const double dX, const double dH) {
	double fxph = f(dX + dH);
	double fxmh = f(dX - dH);

	return (fxph - fxmh) / (2.0 * dH);
}

double HigherDiff(const double (*f)(const double), const double dX, const double dH) {
	double fxpph = f(dX + 2.0 * dH);
	double fxph  = f(dX + dH);
	double fxmh  = f(dX - dH);
	double fxmmh = f(dX - 2.0 * dH);

	return (fxmmh - 8.0 * fxmh + 8.0 * fxph - fxpph) / (12.0 * dH);
}

double ForwardDiff2(const double (*f)(const double), const double dX, const double dH) {
	double fxpph = f(dX + 2.0 * dH);
	double fxph  = f(dX + dH);
	double fx    = f(dX);

	return (fxpph - 2.0 * fxph + fx) / (dH*dH);
}

double BackwardDiff2(const double (*f)(const double), const double dX, const double dH) {
	double fx    = f(dX);
	double fxmh  = f(dX - dH);
	double fxmmh = f(dX - 2.0 * dH);

	return (fxmmh - 2.0 * fxmh + fx) / (dH*dH);

}

double CentralDiff2(const double (*f)(const double), const double dX, const double dH) {
	double fxph = f(dX + dH);
	double fx = f(dX);
	double fxmh = f(dX - dH);

	return (fxph - 2.0 * fx + fxmh) / (dH*dH);
}
