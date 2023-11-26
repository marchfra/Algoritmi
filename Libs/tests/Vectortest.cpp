#include <iostream>
#include <iomanip>
#include <string>

#include "../include/Vector.hpp"

void addVectors(const Vector<double>& v1, const Vector<double>& v2, const Vector<double>& v3);

int main() {
	const int size = 3;
	Vector<double> a("A", size, 0.0);
	Vector<double> b("B", size, 0.0);
	Vector<double> c(size);

	try {
		for (int i = 0; i < size; i++) {
			a[i] = static_cast<double>(i);
			b[i] = static_cast<double>(2 * i);
		}

		std::cout << a << std::endl;
		std::cout << b << std::endl;

		a[5] = 0.0;
	} catch (Vector<double>::Error err) {
		Vector<double>::errorMessage(err);
	} catch (const std::exception& e) {
		std::cout << e.what();
	} catch (...) {
		std::cout << "Uncaught exception" << std::endl;
	}
	return 0;
}
