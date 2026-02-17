#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double f1(double x);
double df1(double x);
double f2(double x);
void bisection(double (*func)(double), double x0, double x1, double tol);
void NewtonRaphson(double (*func)(double), double (*dfunc)(double), double x0, double tol, double dtol);

double f1(double x) {
    return x * x * x - 7. * x * x + 14. * x - 5.;
}

double df1(double x) { // Derivative of f1
    return 3. * x * x - 14. * x + 14.;
}

void bisection(double (*func)(double), double x0, double x1, double tol) {
    double fx0 = func(x0);
    double fx1 = func(x1);
    if (fx0 * fx1 > 0.) {
        printf("Bisection method fails. f(x0) and f(x1) should have opposite signs.\n");
        return;
    }
    double xm, fm;
    int iter = 0;
    while (iter < 1000) {
        xm = 0.5 * (x0 + x1); // Midpoint of the interval
        fm = func(xm);
        if (fx0 * fm < 0.) { // Root is in [x0, xm]
            if (fabs((xm - x0) / xm) < tol) {
            printf("Root found at x = %.15f after %d iterations at relative tolerance %e using the bisection method.\n", xm, iter, tol);
            return; // End the while loop if root is found within tolerance
            }
            x1 = xm;
            fx1 = fm;
        } else { // Root is in [xm, x1]
            if (fabs((xm - x1) / xm) < tol) {
            printf("Root found at x = %.15f after %d iterations at relative tolerance %e using the bisection method.\n", xm, iter, tol);
            return; // End the while loop if root is found within tolerance
            }
            x0 = xm;
            fx0 = fm;
        }
        iter++;
    }
    printf("Bisection method did not converge after 1000 iterations.\n");
}

void NewtonRaphson(double (*func)(double), double (*dfunc)(double), double x0, double tol, double dtol) {
    double epsilon;
    double x1;
    double fx = func(x0);
    double dfx = dfunc(x0);
    int iter = 0;
    while (iter < 10) {
        if (fabs(dfx) < dtol) {
            printf("Derivative is too small. No solution found.\n");
            return; // Avoid division by zero
        }
        x1 = x0 - fx / dfx; // Newton-Raphson update
        epsilon = fabs((x1 - x0) / x1); // Relative error
        // printf("Iteration %d: x = %.15f, f(x) = %.15e, df(x) = %.15e, epsilon = %.15e\n", iter, x1, func(x1), dfunc(x1), epsilon);
        if (epsilon < tol) {
            printf("Root found at x = %.15f after %d iterations at relative tolerance %e using the Newton-Raphson method.\n", x1, iter, tol);
            return; // End the while loop if root is found within tolerance
        }
        x0 = x1;
        fx = func(x0);
        dfx = dfunc(x0);
        iter++;
    }
    printf("Newton-Raphson method did not converge after %d iterations.\n", iter);

}

int main() {
    // Print results for problem 2a
    double tol = 1e-8;
    double dtol = 1e-12;
    bisection(f1, 0., 1., tol);
    NewtonRaphson(f1, df1, 0., tol, dtol);

    return 0;
}