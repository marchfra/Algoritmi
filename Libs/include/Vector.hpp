/**
 * @file       Vector.hpp
 *
 * @brief      This file implements Vector.
 *
 * @author     Francesco Marchisotti.
 *
 * @date       23/11/2023
 */

#pragma once

#include <iostream>
#include <string>
#include <exception>
#include <iomanip>

/**
 * @brief      This class describes a vector.
 *
 * @tparam     T     The tipe of data contained in the vector.
 */
template <class T>
class Vector {
public:
	/**
	 * @brief      This class describes an error.
	 */
	enum class Error {InvalidSize, InvalidIndex, IncompatibleVectors, MemoryAllocation};

	/**
	 * @brief      Constructs a new instance.
	 */
	Vector();
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  size  The size of the vector.
	 */
	Vector(const int& size);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  size          The size.
	 * @param[in]  initialValue  The initial value.
	 */
	Vector(const int& size, T initialValue);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  name  The name of the vector.
	 */
	Vector(const std::string& name);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  name  The name.
	 * @param[in]  size  The size.
	 */
	Vector(const std::string& name, const int& size);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  name          The name.
	 * @param[in]  size          The size.
	 * @param[in]  initialValue  The initial value.
	 */
	Vector(const std::string& name, const int& size, T initialValue);
	/**
	 * @brief      Copy constructor.
	 *
	 * @param[in]  V     The other vector.
	 */
	Vector(const Vector<T>& V);

	~Vector();

	/**
	 * @brief      Sets the size of the vector.
	 *
	 * Used with the default constructor.
	 *
	 * @param[in]  size  The size.
	 */
	void setSize(const int& size);

	/**
	 * @brief      Displays the vector.
	 */
	void display() const;

	// Getter functions
	int getSize() const {
		return m_nRows;
	};
	/**
	 * @brief       Gets the name.
	 *
	 * @param[out]  Name  The name.
	 */
	void getName(std::string& name) const {
		name = m_name;
	};

	// Setter functions
	/**
	 * @brief      Set the value of all elements of a vector.
	 *
	 * @param[in]  value  The value.
	 */
	void set(T value) {
		for (int i = 0; i < m_nRows; i++) m_cells[i] = value;
	}
	/**
	 * @brief      Sets the name.
	 *
	 * @param[in]  Name  The name.
	 */
	void setName(const std::string& name) {
		m_name = name;
	};

	// Overloaded operators
	/**
	 * @brief      Array indexer operator.
	 *
	 * @param[in]  index  The index.
	 *
	 * @return     The element at `index` index.
	 */
	T& operator[](const int& index);
	/**
	 * @brief      Array indexer operator.
	 *
	 * @param[in]  index  The index.
	 *
	 * @return     The element at `index` index.
	 */
	const T& operator[](const int& index) const;
	/**
	 * @brief      Assignment operator.
	 *
	 * Copies all elements of `V` in `this` vector.
	 *
	 * @param[in]  V     The other vector.
	 *
	 * @return     This vector.
	 */
	Vector<T>& operator=(const Vector<T>& V);

	// template <class U>
	friend std::ostream& operator<<(std::ostream& os, const Vector<T>& V);

	static void errorMessage(Vector<T>::Error err) {
		try {
			Vector<T>::errorHandler(err);
		} catch (const std::exception& e) {
			std::cout << e.what();
		}
	}

private:
	int m_nRows;		//< Number of rows in the vector.
	T *m_cells;			//< Address where the vector of type T is stored.
	std::string m_name;	//< The name of the vector.
	/**
	 * @brief      Resource cleanup.
	 */
	void release();
	/**
	 * @brief      Handles errors.
	 *
	 * @param[in]  err   The error.
	 */
	void errorHandler(Vector<T>::Error err) const;
};
