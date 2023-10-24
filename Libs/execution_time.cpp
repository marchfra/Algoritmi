#include <iostream>

#include "execution_time.hpp"

void time_test(void (*f)(), const std::string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	std::cout << f_name << " execution time: " << duration.count() / executions << " ns" << std::endl;
}

void time_test(double (*f)(), const std::string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	std::cout << f_name << " execution time: " << duration.count() / executions << " ns" << std::endl;
}

void time_test(double (*f)(double (*)(double), double, double, int), double(*func)(double), const double a, const double b, int N, const std::string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N);
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	std::cout << f_name << " execution time: " << duration.count() / executions << " ns" << std::endl;
}

void time_test(double (*f)(double (*)(double), double, double, int, int), double(*func)(double), const double a, const double b, int N, int Ng, const std::string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N, Ng);
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	std::cout << f_name << " execution time: " << duration.count() / executions << " ns" << std::endl;
}
