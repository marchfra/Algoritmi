#include <iostream>

using namespace std;

void count_for(int max, bool only_odd);

void count_while(int max, bool only_odd);

int main() {
	int max = 10;
	bool only_odd = true;
	count_while(max, only_odd);

	return 0;
}

void count_for(int max, bool only_odd) {
	for (int i = 0; i < max; i++) {
		if (!(only_odd && i % 2 != 0)) {
			cout << i + 1 << endl;
		}
	}
}

void count_while(int max, bool only_odd) {
	int i = 0;
	while (i < max) {
		if (!(only_odd && i % 2 != 0)) {
			cout << i + 1 << endl;
		}
		i++;
	}
}
