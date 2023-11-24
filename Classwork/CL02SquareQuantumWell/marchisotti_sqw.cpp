// Name: Francesco Marchisotti
// Date: 24/11/2023
// Code output:
// *****************************************************************************
// Bisection, results:
// s = 1; Root = 7.1233753e+00; ntry = 34
// s = 2; Root = 2.8152754e+01; ntry = 34
// s = 3; Root = 6.1433727e+01; ntry = 34

// Secant, results:
// s = 1; Root = 7.1233753e+00; ntry =  6
// s = 2; Root = 2.8152754e+01; ntry =  5
// s = 3; Root = 6.1433727e+01; ntry =  4
// *****************************************************************************

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "../../../Libs/include/root_finder.hpp"

using std::cout;
using std::cin;
using std::endl;

int g_s = 0;
static const double V1 = 250, V2 = 80;

double func(const double& E);

void plot(double (*f)(const double& x));

int main() {
    plot(func);

    double xa, xb;
    const double xtol = 1.0e-8;
    double root;
    int nTry = -1;

    double Vmin;
    if (V1 > V2) Vmin = V2;
    else         Vmin = V1;

    // =============================== Bisection ===============================
    xa = 0.0;   // xa = 0     Domain of sqrt()
    xb = Vmin;  // xb = Vmin  Domain of asin()

    cout << "Bisection, results:" << endl;
    cout.precision(7);
    for (g_s = 1; g_s <= 3; g_s++) {
        bisection(func, xa, xb, xtol, root, nTry);
        cout << "s = " << g_s << "; Root = " << std::scientific << root <<
        "; ntry = " << std::setw(2) << nTry << endl;
    }


    // ================================= Secant ================================
    nTry = -1;

    // The following values have been estimated graphically
    // Arrays with the left and right bounds for the search
    double xL[3] = {0.0, 20.0, 60.0};
    double xR[3] = {10.0, 30.0, 70.0};

    cout << "\nSecant, results:" << endl;
    for (g_s = 1; g_s <= 3; g_s++) {
        secant(func, xL[g_s - 1], xR[g_s - 1], xtol, root, nTry);
        cout << "s = " << g_s << "; Root = " << std::scientific << root <<
        "; ntry = " << std::setw(2) << nTry << endl;
    }

    return 0;
}

double func(const double& E) {
    return static_cast<double>(g_s) * M_PI - asin(sqrt(E / V1))
           - asin(sqrt(E / V2)) - sqrt(E);
}

void plot(double (*f)(const double& x)) {
    std::ofstream out;
    out.open("../data/data.csv");
    if (!out) exit(1);

    const double x0 = 0.0;
    const double xEnd = V2;
    const int nStep = 200;
    const double dx = (xEnd - x0) / nStep;
    double x = x0;

    out << "x,f,s" << endl;
    for (g_s = 1; g_s <= 3; g_s++) {
        for (int i = 0; i < nStep; i++) {
            x = x0 + i * dx;
            out << x << "," << f(x) << "," << g_s << endl;
        }
    }

    out.close();
}


// =====================================================================================================================
// Bisection method
// =====================================================================================================================

// int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
//     int max_ntry = 128;
//     double fa = f(xa);
//     double fb = f(xb);
//     double xm, fm;

//     // Handle fa, fb = 0
//     if (fa == 0.0) {
//         ntry = 0;
//         root = xa;
//         return 0;
//     } else if (fb == 0.0) {
//         ntry = 0;
//         root = xb;
//         return 0;
//     }


//     if (fa * fb < 0) {      // Necessary condition
//         for (int k = 1; k <= max_ntry; k++) {
//             // Function midpoint
//             xm = 0.5 * (xa + xb);
//             fm = f(xm);

//             // Check convergence
//             if (fabs(xb - xa) < xtol || fabs(fm) < ftol || fm == 0.0) {
//                 ntry = k;
//                 root = xm;
//                 return 0;
//             }

//             // Redefine interval
//             if (fm * fa < 0) {
//                 xb = xm;
//                 fb = fm;
//             } else {
//                 xa = xm;
//                 fa = fm;
//             }

//         }

//         ntry = -1;
//         root = nan("");
//         throw std::runtime_error("Maximum number of steps exceeded.");
//     }

//     throw std::runtime_error("The supplied interval does not contain any roots.");
// }

// int bisection(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
//     return bisection(f, xa, xb, xtol, -1.0, root, ntry);
// }

// =====================================================================================================================
// Secant method
// =====================================================================================================================

// int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, const double& ftol, double& root, int& ntry) {
//     int max_ntry = 64;
//     double fa = f(xa);
//     double fb = f(xb);
//     double dx = xb - xa;

//     // Handle fa, fb = 0
//     if (fa == 0.0) {
//         ntry = 0;
//         root = xa;
//         return 0;
//     } else if (fb == 0.0) {
//         ntry = 0;
//         root = xb;
//         return 0;
//     }

//     for (int k = 1; k <= max_ntry; k++) {
//         dx = fb * (xb - xa) / (fb - fa);    // Compute increment

//         // Shift values
//         xa = xb;
//         fa = fb;
//         xb = xb - dx;
//         fb = f(xb);

//         // Check convergence
//         if (fabs(dx) < xtol || fabs(fb) < ftol || fb == 0.0) {
//             ntry = k;
//             root = xb;
//             return 0;
//         }
//     }

//     ntry = -1;
//     root = nan("");
//     throw std::runtime_error("Maximum number of steps exceeded.");
//     return 1;
// }

// int secant(double (*f)(const double& x), double xa, double xb, const double& xtol, double& root, int& ntry) {
//     return secant(f, xa, xb, xtol, -1.0, root, ntry);
// }
