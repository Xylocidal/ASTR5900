#include <stdio.h>
#include <math.h>

double f(double x, double y) {
    return y * y + 1.0;
}

void rk4(double (*f)(double, double), double *x, double *y, double h, int steps) {
    double k1, k2, k3, k4;

    for (int i = 0; i < steps - 1; i++) { // Loop over steps
        k1 = f(x[i], y[i]);
        k2 = f(x[i] + h / 2.0, y[i] + h * k1 / 2.0);
        k3 = f(x[i] + h / 2.0, y[i] + h * k2 / 2.0);
        k4 = f(x[i] + h, y[i] + h * k3);

        y[i + 1] = y[i] + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4); // Update y
        x[i + 1] = x[i] + h; // Update x
    }
}

int main() {
    int steps = 20; // Number of steps
    double h = 0.1; // Step size
    double x[steps]; // Array to store x values
    double y[steps]; // Array to store y values
    x[0] = 0.0; // Initial x value
    y[0] = 0.0; // Initial y value

    rk4(f, x, y, h, steps); // Call the RK4 function

    for(int i = 0; i < steps; i++) {
        printf("x: %.2f, y: %.6f, e: %.6f\n", x[i], y[i], tan(x[i]) - y[i]); // Print results
    }
    return 0;
}