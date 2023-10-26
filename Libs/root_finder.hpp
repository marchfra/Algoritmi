// #ifndef TRUE
// 	#define TRUE 1
// 	#define FALSE 0
// #endif

int bisection(double (*)(double), double, double, const double, const double, double&, int&);

int false_position(double (*)(double), double, double, const double, const double, double&, int&);

int secant(double (*)(double), double, double, const double, const double, double&, int&);

int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&, int&);
