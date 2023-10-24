#include <iostream>

using namespace std;

void intrand();

void doublerand();

void time_test(void (*)(), const string, const long int = 1e6);

int main() {
	long int n = 1e8;

	srand(time(NULL));
	srand48(time(NULL));

	time_test(intrand, 	"intrand()", n);
	time_test(doublerand, "doublerand()", n);

	return 0;
}

void time_test(void (*f)(), const string f_name, const long int executions) {
	auto start = std::chrono::high_resolution_clock::now();
	for (long int i = 0; i < executions; i++) {
		f();
	}
	auto stop = std::chrono::high_resolution_clock::now();	
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

	cout << f_name << " execution time: " << duration.count() / executions << " ns" << endl;
}


void intrand() {
	bool x = (rand() % 2 == 0);
}

void doublerand() {
	bool x = (drand48() >= 0.5);
}
