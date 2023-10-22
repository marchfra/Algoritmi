#include <iostream>

using namespace std;

void add_one(int&);

int main() {
	int a = 5;

	add_one(a);

	cout << a << endl;

	return 0;
}

void add_one(int& a) {
	a += 1;
}
