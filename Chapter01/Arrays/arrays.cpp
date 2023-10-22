#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

#define NSIZE 1000

void stat(int[], int, double&, double&, double&);

int main() {
	srand(time(NULL));

	int size = NSIZE, maxVal = 100;
	int a[size];
	double mu, s2, s;

	for (int i = 0; i < size; i++) {
		int randNum = rand();
		a[i] = randNum % maxVal;
	}

	stat(a, size, mu, s2, s);

	cout << mu << " " << s2 << " " << s << endl;

	return 0;
}

void stat(int arr[], int size, double& mu, double& s2, double& s) {
	int sum = 0;
	for (int i = 0; i < size; i++) {
		sum += arr[i];
	}
	mu = (double)sum / size;
	
	s2 = 0.0;
	for (int i = 0; i < size; i++) {
		s2 += (arr[i] - mu) * (arr[i] - mu);
	}

	s = sqrt(s2);
}
