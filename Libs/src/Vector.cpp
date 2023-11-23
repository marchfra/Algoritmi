#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "Vector.hpp"


// ================================ Constructors ===============================
template <class T>
Vector<T>::Vector() {
	nRows = 0;
	cells = nullptr;
}

template <class T>
Vector<T>::Vector(int size) {
	setSize(size);
	cells = nullptr;
}

template <class T>
Vector<T>::Vector(int size, T initialValue) {
	setSize(size);
	set(initialValue);
}

template <class T>
Vector<T>::Vector(const std::string& name) {`
	setName(name);
	cells = nullptr;
}

template <class T>
Vector<T>::Vector(const std::string& name, int size) {
	setName(name);
	setSize(size);
	cells = nullptr;
}

template <class T>
Vector<T>::Vector(const std::string& name, int size, T initialValue) {
	setName(name);
	setSize(size);
	set(initialValue);
}

template <class T>
Vector<T>::Vector(const Vector<T>& V) {
	setName(V.name);
	setSize(V.getSize());
	for (int i = 0; i < nRows; i++) {
		cells[i] = V[i];
	}
}

// ================================ Destructors ================================
template <class T>
Vector<T>::~Vector() {
	release();
}

// ============================== Helper functions =============================
template <class T>
void Vector<T>::setSize(int size) {
	if (size <= 0) {
		errorHandler(Vector<T>::error::InvalidSize);
	}
	release();
	try {
		cells = new T[size];
	} catch (std::bad_alloc& e) {
		errorHandler(Vector<T>::error::MemoryAllocation);
	}

	nRows = size;
}

template <class T>
Vector<T>::release() {
	if (cells != nullptr) {
		delete[] cells;
		cells = nullptr;
		nRows = 0;
	}
}

// ============================ Overloaded operators ===========================
template <class T>
T& Vector<T>::operator[](int index) {
	if (index < 0 || index >= nRows) {
		errorHandler(Vector<T>::error::InvalidIndex);
		return cells[0];
	} else {
		return cells[index];
	}
}

template <class T>
const T& Vector<T>::operator[](int index) const {
	if (index < 0 || index >= nRows) {
		errorHandler(Vector<T>::error::InvalidIndex);
		return cells[0];
	} else {
		return cells[index];
	}
}

template <class T>
T& Vector<T>::operator=(const Vector<T>& V) {
	if (this != &V) {	// Protects against self-assignment
		if (nRows != V.getSize()) {	// Compatible vectors?
			errorHandler(Vector<T>::error::IncompatibleVectors);
			return *this;
		}
		for (int i = 0; i < nRows; i++) {
			cells[i] = V[i];
		}
	}

	return *this;
}

// =============================== Error handler ===============================
template <class T>
void Vector<T>::errorHandler(Error err) const {
	std::string errors[] = {"Invalid vector size to create vector.",
							"Invalid index to access vector element.",
							"Vector operation cannot take place due to incompatible vectors.",
							"Vector cannot be created. Memory allocation error."};

	throw std::exception(errors[err]);
}
