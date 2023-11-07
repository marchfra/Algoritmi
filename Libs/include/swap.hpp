#pragma once

// void swap(int&, int&);

// void swap(double&, double&);

template <class T>
void swap(T& a, T& b) {
	T temp = b;
	b = a;
	a = temp;
}
