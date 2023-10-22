#include <iostream>

using namespace std;

int quotient(int, int, int&, int&);

int main() {
	int a = 13, b = 6;
	int quot, rem;

	int flag = quotient(a, b, quot, rem);

	if (flag == 1) {
		cout << "Division by 0!" << endl;
	} else if (flag == 0) {
		cout << a << " / " << b << " = " << quot << " remainder " << rem << endl;	
	}


	return 0;
}

int quotient(int a, int b, int& quot, int& rem) {
	if (b == 0) return 1;	// Failure
	quot = a / b;
	rem = a % b;
	return 0;
}
