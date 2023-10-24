#include <iostream>
#include "swap.h"

using namespace std;

int main() {
	int a = 1, b = 13;

	cout << "a = " << a <<"\tb = " << b << endl;
	swap(a, b);
	cout << "a = " << a <<"\tb = " << b << endl;
}
