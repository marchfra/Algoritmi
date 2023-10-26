#include <iostream>

using namespace std;

#include "debug.hpp"
#include "execution_time.hpp"

void time_test(void (*f)(), const string f_name, const long int executions) {
	auto start = chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = chrono::high_resolution_clock::now();	
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(), const string f_name, const long int executions) {
	auto start = chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = chrono::high_resolution_clock::now();	
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(double (*)(double), double, double, int), double(*func)(double), const double a, const double b, int N, const string f_name, const long int executions) {
	auto start = chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N);
	}
	auto stop = chrono::high_resolution_clock::now();	
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}

void time_test(double (*f)(double (*)(double), double, double, int, int), double(*func)(double), const double a, const double b, int N, int Ng, const string f_name, const long int executions) {
	auto start = chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f(func, a, b, N, Ng);
	}
	auto stop = chrono::high_resolution_clock::now();	
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}
