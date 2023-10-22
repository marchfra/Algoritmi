#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

double heron(double, double, double = 1.0e-16);

int main() {

	cout << setiosflags(ios::scientific);

	double S, guess;

	cout << "Enter a non-negative real number: ";
	do {
		cin >> S;
		if (S < 0) cout << "You entered a negative number. Try again!\n\
							Enter a non-negative real number: ";
	} while (S < 0);

	cout << "Enter yor guess: ";
	do {
		cin >> guess;
		if (guess <= 0) cout << "You entered a non-positive number. Try again!\n\
								 Enter a non-negative real number: ";
	} while (guess <= 0);

	cout << "=========================" << endl;

	double root = heron(S, guess);

	cout << endl;
	cout << "The estimated sqrt of " << S << " is: " << root << endl;
	cout << "The exact value is:                    " << sqrt(S) << endl;

	return 0;
}

double heron(double S, double x_n, double tol) {
	double x_N, err;
	int i = 0;

	do {
		x_N = 0.5 * (x_n + S / x_n);
		err = fabs(x_N - x_n);

		cout << "Iteration #" << setw(2) << setfill(' ') << i + 1 << "; ";
		cout << "x = " << x_N << "; err = " << err << endl;

		x_n = x_N;
		i++;
	} while (err > 0);

	return x_N;
}
