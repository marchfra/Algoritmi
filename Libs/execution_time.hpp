void time_test(void (*)(), const std::string, const long int = 1e6);

void time_test(double (*)(), const std::string, const long int = 1e6);

// Time test for Rectangular, Midpoint, Trapezoid, Simpson quad
void time_test(double (*)(double (*)(double), double, double, int), double(*)(double), const double, const double, int, const std::string, const long int = 1e6);

// Time test for Gauss quad
void time_test(double (*)(double (*)(double), double, double, int, int), double(*)(double), const double, const double, int, int, const std::string, const long int = 1e6);
