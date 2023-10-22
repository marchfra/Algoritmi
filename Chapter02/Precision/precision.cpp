#include <iostream>
#include <iomanip>

using namespace std;

int main() {
	float oneF = 1.0, epsF = 1.0;
	double oneD = 1.0, epsD = 1.0;
	long double oneL = 1.0, epsL = 1.0;
	int iterF = 0, iterD = 0, iterL = 0;

	while (oneF + epsF != oneF) {	// while (1.0 + eps != 1.0) yields the wrong result,
									// because 1.0 is (almost) always interpreted as a double
		iterF++;
		epsF /= 10;
	}
	cout << "Float precision:        " << iterF << "; eps = " << epsF << endl;

	while (oneD + epsD != oneD) {
		iterD++;
		epsD /= 10;
	}
	cout << "Double precision:      " << iterD << "; eps = " << epsD << endl;

	while (oneL + epsL != oneL) {
		iterL++;
		epsL /= 10;
	}
	cout << "Long double precision: " << iterL << "; eps = " << epsL << endl << endl;


	const float PI_f = 3.141592653589793238462643383279502884197;

	cout << "float precision:" << endl;
	cout << "         ↓" << endl;
	cout << "3.141592653589793238462643383279502884197" << endl;
	// cout << setprecision(40) << PI_f << endl;
	printf("%.39f\n", PI_f);
	cout << "        ↑" << endl;
	cout << "Actual precision:" << endl << endl;

	const double PI_d = 3.141592653589793238462643383279502884197;

	cout << "double precision:" << endl;
	cout << "                 ↓" << endl;
	cout << "3.141592653589793238462643383279502884197" << endl;
	// cout << setprecision(40) << PI_d << endl;
	printf("%.39lf\n", PI_d);
	cout << "                 ↑" << endl;
	cout << "Actual precision:" << endl << endl;

	const long double PI_ld = 3.141592653589793238462643383279502884197;

	cout << "long double precision:" << endl;
	cout << "                     ↓" << endl;
	cout << "3.141592653589793238462643383279502884197" << endl;
	// cout << setprecision(40) << PI_ld << endl;
	printf("%.39Lf\n", PI_ld);
	cout << "                 ↑" << endl;
	cout << "Actual precision:" << endl << endl;

	// cout << "Size of bool: " << sizeof(bool) << "B" << endl;
	// cout << "Size of unsigned char: " << sizeof(unsigned char) << "B" << endl;
	// cout << "Size of char: " << sizeof(char) << "B" << endl;
	// cout << "Size of wchar_t: " << sizeof(wchar_t) << "B" << endl;
	// cout << "Size of short: " << sizeof(short) << "B" << endl;
	// cout << "Size of unsigned short: " << sizeof(unsigned short) << "B" << endl;
	// cout << "Size of int: " << sizeof(int) << "B" << endl;
	// cout << "Size of unsigned int: " << sizeof(unsigned int) << "B" << endl;
	// cout << "Size of long: " << sizeof(long) << "B" << endl;
	// cout << "Size of unsigned long: " << sizeof(unsigned long) << "B" << endl;
	// cout << "Size of float: " << sizeof(float) << "B" << endl;
	// cout << "Size of long long: " << sizeof(long long) << "B" << endl;
	// cout << "Size of unsigned long long: " << sizeof(unsigned long long) << "B" << endl;
	// cout << "Size of double: " << sizeof(double) << "B" << endl;
	// cout << "Size of long double: " << sizeof(long double) << "B" << endl;

	return 0;
}
