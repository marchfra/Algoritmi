#include <iostream>

using namespace std;

int array_size(double arr[]);

int main() {
	double arr[15];

	int nSize = sizeof(arr) / sizeof(arr[0]);

	cout << "main()" << endl;
	cout << "sizeof(arr)    = " << sizeof(arr) << endl;
	cout << "sizeof(arr[0]) = " << sizeof(arr[0]) << endl;
	cout << "nSize from main: " << nSize << endl << endl;

	array_size(arr);

	return 0;
}

int array_size(double arr[]) {
	int nSize = sizeof(*arr) / sizeof(arr[0]);
	cout << "array_size()" << endl;
	cout << "sizeof(arr)    = " << sizeof(arr) << endl;
	cout << "sizeof(*arr)   = " << sizeof(*arr) << endl;
	cout << "sizeof(arr[0]) = " << sizeof(arr[0]) << endl;
	cout << "nSize from func: " << nSize << endl;
	return nSize;
}
