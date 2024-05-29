/**
 * @file    exception.hpp
 *
 * @author  Francesco Marchisotti
 *
 * @brief   Implements the class exception
 *
 * @date    2024-05-08
 */
#pragma once

#include <iostream>

class exception : public std::exception {
public:
  /**
   * @brief Constructor (C++ STL strings).
   *
   * @param message The error message.
   */
  explicit exception(const std::string& message)
    : msg_(message) {}

  /**
   * @brief Destructor.
   *
   * Virtual to allow for subclassing.
   */
  virtual ~exception() noexcept {}

  /**
   * @brief Returns a pointer to the (constant) error description.
   *
   * @return A pointer to a const char*. The underlying memory is in posession
   *         of the exception object. Callers must not attempt to free the
   *         memory.
   */
  virtual const char* what() const noexcept { return msg_.c_str(); }

protected:
  std::string msg_;  //!< Error message
};
