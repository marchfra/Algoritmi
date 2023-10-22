#include <iostream>

using namespace std;

int sum(int, int);

int main() {
	int a = 15, b = 18;

	cout << a << " + " << b << " = " << sum(a, b) << endl;

	return 0;
}

int sum(int a, int b) {
	return a + b;
}
