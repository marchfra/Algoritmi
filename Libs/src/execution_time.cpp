#include "../include/execution_time.hpp"

#include "../include/debug.hpp"

double timeTest(void (*f)(), const std::string &f_name,
                const long int executions, bool print) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) f();
	auto stop = std::chrono::high_resolution_clock::now();

	auto duration =
		std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	double exec_time = duration.count() / executions;
	if (print)
		std::cout << f_name << " execution time: " << exec_time << " ns"
				  << std::endl;

	return exec_time;
}

double timeTest(void (*f)(), const long int executions, bool print) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) f();
	auto stop = std::chrono::high_resolution_clock::now();

	auto duration =
		std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	double exec_time = duration.count() / executions;
	if (print)
		std::cout << "execution time: " << exec_time << " ns" << std::endl;

	return exec_time;
}
