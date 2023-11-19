/**
 * @file execution_time.hpp
 *
 * @brief      Implementetion of execution time utility functions.
 * 
 * @author     Francesco Marchisotti
 * 
 * @date       19/11/2023
 */
#pragma once

#include <iostream>

/**
 * @brief      Measure execution time of a function.
 *
 * Measure average execution time of a function over a number of
 * executions.
 *
 * @param[in]  f           The function to be measured.
 * @param[in]  fName       The name of the function.
 * @param[in]  executions  The number of executions.
 * @param[in]  print       Print the execution time in addition to returning it.
 *
 * @return     Average execution time.
 */
double timeTest(void (*f)(), const std::string& fName, const long int executions = 1e6, bool print = false);

/**
 * @overload
 *
 * @brief      Measure execution time of a function.
 *
 * Measure average execution time of a function over a number of
 * executions.
 *
 * @param[in]  f           The function to be measured.
 * @param[in]  executions  The number of executions.
 * @param[in]  print       Print the execution time in addition to returning it.
 *
 * @return     Average execution time.
 */
double timeTest(void (*f)(), const long int executions = 1e6, bool print = false);
