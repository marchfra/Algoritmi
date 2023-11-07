#include <vector>

LegendrePolynomial::LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations)
    : mLowerBound(lowerBound), mUpperBound(upperBound), mNumberOfIterations(numberOfIterations), mWeight(numberOfIterations+1), mRoot(numberOfIterations+1) {
    calculateWeightAndRoot();
}

const std::vector<double>& LegendrePolynomial::getWeight() const {
    return mWeight;
}

const std::vector<double>& LegendrePolynomial::getRoot() const {
    return mRoot;
}

void LegendrePolynomial::calculateWeightAndRoot() {
    for (int step = 0; step <= mNumberOfIterations; step++) {
        double root = cos(M_PI * (step - 0.25) / (mNumberOfIterations + 0.5));
        Result result = calculatePolynomialValueAndDerivative(root);

        double newtonRaphsonRatio;
        do {
            newtonRaphsonRatio = result.value / result.derivative;
            root -= newtonRaphsonRatio;
            result = calculatePolynomialValueAndDerivative(root);
        } while (fabs(newtonRaphsonRatio) > EPSILON);

        mRoot[step] = root;
        mWeight[step] = 2.0 / ((1 - root*root) * result.derivative * result.derivative);
    }
}

Result LegendrePolynomial::calculatePolynomialValueAndDerivative(double x) {
    Result result(x, 0);

    double value_minus_1 = 1;
    const double f = 1 / (x*x - 1);
    for (int step = 2; step <= mNumberOfIterations; step++) {
        const double value = ((2 * step - 1) * x * result.value - (step - 1) * value_minus_1) / step;
        result.derivative = step * f * (x * value - result.value);

        value_minus_1 = result.value;
        result.value = value;
    }

    return result;
}

const double LegendrePolynomial::EPSILON = 1e-15;

// double gaussLegendreIntegral(double a, double b, int n, const std::function<double (double)>& f) {
//     const LegendrePolynomial legendrePolynomial(a, b, n);
//     const std::vector<double>& weight = legendrePolynomial.getWeight();
//     const std::vector<double>& root = legendrePolynomial.getRoot();

//     const double width = 0.5 * (b - a);
//     const double mean  = 0.5 * (a + b);

//     double gaussLegendre = 0;
//     for(int step = 1; step <= n; step++) {
//         gaussLegendre += weight[step] * f(width * root[step] + mean);
//     }

//     return gaussLegendre * width;
// }
