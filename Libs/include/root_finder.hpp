#pragma once

#include <iostream>
#include <iomanip>

int find_roots(double (*)(double), double (*)(double), const double, const double, const double, double[], int&, const int = 128, const std::string = "newton");
int find_roots(double (*)(double), const double, const double, const double, double[], int&, const int = 128, const std::string = "bisection");

void bracket(double (*)(double), const double, const double, double[], double[], const int, int&);

int bisection(double (*)(double), double, double, const double, const double, double&, int&);
int bisection(double (*)(double), double, double, const double, double&);
int bisection(double (*)(double), double, double, const double, double&, int&);
int bisection(double (*)(double), double, double, const double, const double, double&);

int false_position(double (*)(double), double, double, const double, const double, double&, int&);
int false_position(double (*)(double), double, double, const double, double&);
int false_position(double (*)(double), double, double, const double, double&, int&);
int false_position(double (*)(double), double, double, const double, const double, double&);

int secant(double (*)(double), double, double, const double, const double, double&, int&);
int secant(double (*)(double), double, double, const double, double&);
int secant(double (*)(double), double, double, const double, double&, int&);
int secant(double (*)(double), double, double, const double, const double, double&);

int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&, int&);
int newton(double (*)(double), double (*)(double), double, double, const double, double&);
int newton(double (*)(double), double (*)(double), double, double, const double, double&, int&);
int newton(double (*)(double), double (*)(double), double, double, const double, const double, double&);
