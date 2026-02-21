#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
    double omega_0;
    double a;
    double beta;
} ParamsF2;

double f1_ctx(double x, void *ctx);
double df1_ctx(double x, void *ctx);

double f2_ctx(double omega, void *ctx);
double df2_ctx(double omega, void *ctx);

void bisection_ctx(double (*func)(double, void *), void *ctx,
                   double x0, double x1, double tol, bool info_toggle);

void NewtonRaphson_ctx(double (*func)(double, void *), double (*dfunc)(double, void *),
                       void *ctx, double x0, double tol, double dtol, bool info_toggle);

double f1_ctx(double x, void *ctx) {
    (void)ctx;
    return x * x * x - 7. * x * x + 14. * x - 5.;
}

double df1_ctx(double x, void *ctx) {
    (void)ctx;
    return 3. * x * x - 14. * x + 14.;
}

double f2_ctx(double omega, void *ctx) {
    ParamsF2 *p = (ParamsF2 *)ctx;
    double omega_0 = p->omega_0;
    double a       = p->a;
    double beta    = p->beta;

    return (2. * (omega - omega_0)) / (a * a + (omega - omega_0) * (omega - omega_0))
         + beta / (exp(beta * omega) - 1.);
}

double df2_ctx(double omega, void *ctx) {
    ParamsF2 *p = (ParamsF2 *)ctx;
    double omega_0 = p->omega_0;
    double a       = p->a;
    double beta    = p->beta;

    double d = omega - omega_0;
    double denom = a * a + d * d;

    // term from beta/(exp(beta*omega)-1)
    double e = exp(beta * omega);
    double term1 = -(beta * beta * e) / ((e - 1.) * (e - 1.));

    // derivative of 2d/(a^2+d^2)
    double term2 = 2. / denom;
    double term3 = -4. * d * d / (denom * denom);

    return term1 + term2 + term3;
}

void bisection_ctx(double (*func)(double, void *), void *ctx,
                   double x0, double x1, double tol, bool info_toggle) {

    double fx0 = func(x0, ctx);
    double fx1 = func(x1, ctx);

    if (fx0 * fx1 > 0.) {
        printf("Bisection method fails. f(x0) and f(x1) should have opposite signs.\n");
        return;
    }

    int iter = 1;
    while (iter < 1000) {
        double xm = 0.5 * (x0 + x1);
        double fm = func(xm, ctx);

        double epsilon = fabs(x1 - x0) / fmax(1.0, fabs(xm)); // interval-based criterion

        if (info_toggle) {
            printf("Iteration %d: x = %.15f\n", iter, xm);
        }
        if (epsilon < tol) {
            printf("Root found at x = %.15f after %d iterations\n"
                   "at relative tolerance %e using the bisection method.\n",
                   xm, iter, tol);
            return;
        }

        if (fx0 * fm < 0.) { x1 = xm; fx1 = fm; }
        else               { x0 = xm; fx0 = fm; }

        iter++;
    }

    printf("Bisection method did not converge after 1000 iterations.\n");
}

void NewtonRaphson_ctx(double (*func)(double, void *), double (*dfunc)(double, void *),
                       void *ctx, double x0, double tol, double dtol, bool info_toggle) {
    double x1, epsilon;
    int iter = 1;

    double fx = func(x0, ctx);
    double dfx = dfunc(x0, ctx);

    while (iter < 1000) {
        if (fabs(dfx) < dtol) {
            printf("Derivative is too small. Loop terminated.\n");
            return;
        }

        x1 = x0 - fx / dfx;

        // safer relative error in case x1~0
        epsilon = fabs(x1 - x0) / fmax(1.0, fabs(x1));

        if (info_toggle) {
            printf("Iteration %d: x = %.15f\n", iter, x1);
        }
        if (epsilon < tol) {
            printf("Root found at x = %.15f after %d iterations\n"
                   "at relative tolerance %e using the Newton-Raphson method.\n",
                   x1, iter, tol);
            return;
        }

        x0 = x1;
        fx = func(x0, ctx);
        dfx = dfunc(x0, ctx);
        iter++;
    }

    printf("Newton-Raphson method did not converge after %d iterations.\n", iter);
}

int main() {
    // Print results for problem 2a
    printf("Problem 2a:\n");
    double tol = 1e-8;
    double dtol = 1e-12;
    bisection_ctx(f1_ctx, NULL, 0., 1., tol, true);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 0., tol, dtol, true);

    // Print results for problem 2b
    printf("\nProblem 2b:\n");
    bisection_ctx(f1_ctx, NULL, 0., 1., tol, false);
    bisection_ctx(f1_ctx, NULL, 0.3, 0.75, tol, false);
    bisection_ctx(f1_ctx, NULL, 0.1, 2., tol, false);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 0., tol, dtol, false);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 0.5, tol, dtol, false);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 1., tol, dtol, false);

    // Print results for problem 2c
    printf("\nProblem 2c:\n");
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 1.5, tol, dtol, false);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 2., tol, dtol, false);
    NewtonRaphson_ctx(f1_ctx, df1_ctx, NULL, 3.21525, tol, dtol, false);

    // Print results for problem 3
    printf("\nProblem 3:\n");
    ParamsF2 p1 = {.omega_0 = 50.0, .a = 15.0, .beta = 0.02};
    ParamsF2 p2 = {.omega_0 = 50.0, .a = 15.0, .beta = 0.05};
    ParamsF2 p3 = {.omega_0 = 50.0, .a = 15.0, .beta = 0.1};

    NewtonRaphson_ctx(f2_ctx, df2_ctx, &p1, 50.0, tol, dtol, false); // initial guess at resonance
    NewtonRaphson_ctx(f2_ctx, df2_ctx, &p2, 50.0, tol, dtol, false);
    NewtonRaphson_ctx(f2_ctx, df2_ctx, &p3, 50.0, tol, dtol, false);

    return 0;
}