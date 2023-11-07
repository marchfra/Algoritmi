#pragma once

void EulerStep(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq);

void RK2Step(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq, bool rule = false);

void RK4Step(double t, double Y[], void (*RHSFunc)(double t, double Y[], double R[]), double dt, int neq);
