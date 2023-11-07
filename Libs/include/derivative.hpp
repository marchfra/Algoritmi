#pragma once

double ForwardDiff(const double (*f)(const double), const double dX, const double dH);

double BackwardDiff(const double (*f)(const double), const double dX, const double dH);

double CentralDiff(const double (*f)(const double), const double dX, const double dH);

double HigherDiff(const double (*f)(const double), const double dX, const double dH);

double ForwardDiff2(const double (*f)(const double), const double dX, const double dH);

double BackwardDiff2(const double (*f)(const double), const double dX, const double dH);

double CentralDiff2(const double (*f)(const double), const double dX, const double dH);
