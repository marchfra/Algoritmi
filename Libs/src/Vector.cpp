// #include <iostream>
// #include <iomanip>
// #include <stdexcept>

#include "../include/Vector.hpp"


// ================================ Constructors ===============================
template <class T>
Vector<T>::Vector() {
	m_nRows = 0;
	m_cells = nullptr;
}

template <class T>
Vector<T>::Vector(const int& size) {
	setSize(size);
	m_cells = nullptr;
}

template <class T>
Vector<T>::Vector(const int& size, T initialValue) {
	setSize(size);
	set(initialValue);
}

template <class T>
Vector<T>::Vector(const std::string& name) {
	setName(name);
	m_cells = nullptr;
}

template <class T>
Vector<T>::Vector(const std::string& name, const int& size) {
	setName(name);
	setSize(size);
	m_cells = nullptr;
}

template <class T>
Vector<T>::Vector(const std::string& name, const int& size, T initialValue) {
	setName(name);
	setSize(size);
	set(initialValue);
}

template <class T>
Vector<T>::Vector(const Vector<T>& V) {
	setName(V.name);
	setSize(V.getSize());
	for (int i = 0; i < m_nRows; i++) {
		m_cells[i] = V[i];
	}
}

// ================================ Destructors ================================
template <class T>
Vector<T>::~Vector() {
	release();
}

// ============================== Helper functions =============================
template <class T>
void Vector<T>::setSize(const int& size) {
	if (size <= 0) {
		errorHandler(Vector<T>::Error::InvalidSize);
	}
	release();
	try {
		m_cells = new T[size];
	} catch (std::bad_alloc& e) {
		errorHandler(Vector<T>::Error::MemoryAllocation);
	}

	m_nRows = size;
}

template <class T>
void Vector<T>::release() {
	if (m_cells != nullptr) {
		delete[] m_cells;
		m_cells = nullptr;
		m_nRows = 0;
	}
}

// ============================ Overloaded operators ===========================
template <class T>
T& Vector<T>::operator[](const int& index) {
	if (index < 0 || index >= m_nRows) {
		errorHandler(Vector<T>::Error::InvalidIndex);
		return m_cells[0];
	} else {
		return m_cells[index];
	}
}

template <class T>
const T& Vector<T>::operator[](const int& index) const {
	if (index < 0 || index >= m_nRows) {
		errorHandler(Vector<T>::Error::InvalidIndex);
		return m_cells[0];
	} else {
		return m_cells[index];
	}
}

template <class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& V) {
	if (this != &V) {	// Protects against self-assignment
		if (m_nRows != V.getSize()) {	// Compatible vectors?
			errorHandler(Vector<T>::Error::IncompatibleVectors);
			return *this;
		}
		for (int i = 0; i < m_nRows; i++) {
			m_cells[i] = V[i];
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
