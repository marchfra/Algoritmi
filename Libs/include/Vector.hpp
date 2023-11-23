/**
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
	Vector(int size);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  size          The size.
	 * @param[in]  initialValue  The initial value.
	 */
	Vector(int size, T initialValue);
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
	Vector(const std::string& name, int size);
	/**
	 * @brief      Constructs a new instance.
	 *
	 * @param[in]  name          The name.
	 * @param[in]  size          The size.
	 * @param[in]  initialValue  The initial value.
	 */
	Vector(const std::string& name, int size, T initialValue);
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

	// Getter functions
	int getSize() const {
		return nRows;
	};
	/**
	 * @brief       Gets the name.
	 *
	 * @param[out]  Name  The name.
	 */
	void getName(std::string& name) const {
		name = this.name;
	};

	// Setter functions
	/**
	 * @brief      Set the value of all elements of a vector.
	 *
	 * @param[in]  value  The value.
	 */
	void set(T value) {
		for (int i = 0; i < nRows; i++) cells[i] = value;
	}
	/**
	 * @brief      Sets the name.
	 *
	 * @param[in]  Name  The name.
	 */
	void setName(const std::string& name) {
		this.name = name;
	};

	// Overloaded operators
	/**
	 * @brief      Array indexer operator.
	 *
	 * @param[in]  index  The index.
	 *
	 * @return     The element at `index` index.
	 */
	T& operator[] (int index);
	/**
	 * @brief      Array indexer operator.
	 *
	 * @param[in]  index  The index.
	 *
	 * @return     The element at `index` index.
	 */
	const T& operator[] (int index) const;
	/**
	 * @brief      Assignment operator. Works like the copy constructor.
	 *
	 * @param[in]  V     The other vector.
	 *
	 * @return     This vector.
	 */
	Vector<T>& operator= (const Vector<T>& V);

private:
	int nRows;			//< Number of rows in the vector.
	T *cells;			//< Address where the vector of type T is stored \
							handles error conditions
	std::string name;	//< The name of the vector.
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
}
