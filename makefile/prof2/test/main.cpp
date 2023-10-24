#include <iostream>
#include "../Libs/swap.hpp"

using namespace std;

int main() {
	cout << "Hello World!" << endl;

	int a = 1, b = 32;

	cout << "a = " << a << "\tb = " << b << endl;

	cout << "Swap..." << endl;

	swap(a, b);

	cout << "a = " << a << "\tb = " << b << endl;

	return 0;
}
