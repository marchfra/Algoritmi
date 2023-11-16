/**
 * @file execution_time.hpp
 * 
 * @brief Implementetion of execution time utility functions.
 */
#pragma once

#include <iostream>

/**
 * @brief Measure execution time of a function.
 * 
 * Measure average execution time of a function over a number of executions.
 * 
 * @param[in] f The function to be measured.
 * @param[in] f_name The name of the function.
 * @param[in] executions The number of executions.
 * @param[in] print Print the execution time in addition to returning it.
 * 
 * @returns Average execution time.
 */
double time_test(void (*)(), const std::string&, const long int = 1e6, bool = false);

/**
 * @overload
 * 
 * @brief Measure execution time of a function.
 * 
 * Measure average execution time of a function over a number of executions.
 * 
 * @param[in] f The function to be measured.
 * @param[in] executions The number of executions.
 * @param[in] print Print the execution time in addition to returning it.
 * 
 * @returns Average execution time.
 */
double time_test(void (*)(), const long int = 1e6, bool = false);
