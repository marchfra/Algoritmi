#pragma once

#include <vector>

class LegendrePolynomial {
public:
    LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations);

    const std::vector<double>& getWeight() const;

    const std::vector<double>& getRoot() const;

private:
    const static double EPSILON;

    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    void calculateWeightAndRoot();

    Result calculatePolynomialValueAndDerivative(double x);

    const double mLowerBound;
    const double mUpperBound;
    const int mNumberOfIterations;
    std::vector<double> mWeight;
    std::vector<double> mRoot;
};

// double gaussLegendreIntegral(double a, double b, int n, const std::function<double (double)>& f);
