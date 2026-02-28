#include <stdio.h>
#include <math.h>

double f(double x, double y);
void rk4(double (*f)(double, double), double *x, double *y, double h, int steps);
void euler(double (*f)(double, double), double *x, double *y, double h, int steps);

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

void euler(double (*f)(double, double), double *x, double *y, double h, int steps) {
    for (int i = 0; i < steps - 1; i++) { // Loop over steps
        y[i + 1] = y[i] + h * f(x[i], y[i]); // Update y
        x[i + 1] = x[i] + h; // Update x
    }
}

int main() {
    int steps = 21; // Number of steps
    int steps2 = 41; // More steps for smaller h
    int steps3 = 201; // Even more steps for smaller h
    int steps4 = 401; // Even more steps for smaller h
    double h = 0.1; // Step size
    double h2 = 0.05; // Smaller step size for comparison
    double h3 = 0.01; // Even smaller step size for comparison
    double h4 = 0.005; // Even smaller step size for comparison
    double x[steps]; // Array to store x values
    double y[steps]; // Array to store y values
    double x2[steps2]; // Array for smaller step size
    double y2[steps2]; // Array for smaller step size
    double x3[steps3]; // Array for even smaller step size
    double y3[steps3]; // Array for even smaller step size
    double x4[steps4]; // Array for even smaller step size
    double y4[steps4]; // Array for even smaller step size
    x[0] = 0.0; // Initial x value
    y[0] = 0.0; // Initial y value

    // Perform Euler's method with h = 0.1
    euler(f, x, y, h, steps);

    FILE *fp = fopen("euler_results_h_0.1.txt", "w"); // Print results to file
    for(int i = 0; i < steps; i++) {
        fprintf(fp, "%.15f\t%.15f\t%.15f\n", x[i], y[i], tan(x[i]) - y[i]);
    }
    fclose(fp);

    x2[0] = 0.0; // Reset initial x value
    y2[0] = 0.0; // Reset initial y value

    // Perform Euler's method with h = 0.05
    euler(f, x2, y2, h2, steps2);

    FILE *fp2 = fopen("euler_results_h_0.05.txt", "w");
    for(int i = 0; i < steps2; i++) {
        fprintf(fp2, "%.15f\t%.15f\t%.15f\n", x2[i], y2[i], tan(x2[i]) - y2[i]);
    }
    fclose(fp2);

    x3[0] = 0.0; // Reset initial x value
    y3[0] = 0.0; // Reset initial y value

    // Perform Euler's method with h = 0.01
    euler(f, x3, y3, h3, steps3);

    FILE *fp3 = fopen("euler_results_h_0.01.txt", "w");
    for(int i = 0; i < steps3; i++) {
        fprintf(fp3, "%.15f\t%.15f\t%.15f\n", x3[i], y3[i], tan(x3[i]) - y3[i]);
    }
    fclose(fp3);
    
    x4[0] = 0.0; // Reset initial x value
    y4[0] = 0.0; // Reset initial y value

    // Perform Euler's method with h = 0.005
    euler(f, x4, y4, h4, steps4);

    FILE *fp4 = fopen("euler_results_h_0.005.txt", "w");
    for(int i = 0; i < steps4; i++) {
        fprintf(fp4, "%.15f\t%.15f\t%.15f\n", x4[i], y4[i], tan(x4[i]) - y4[i]);
    }
    fclose(fp4);
    
    x[0] = 0.0; // Reset initial x value
    y[0] = 0.0; // Reset initial y value

    // Perform RK4 method with h = 0.1
    rk4(f, x, y, h, steps);

    FILE *fp5 = fopen("rk4_results_h_0.1.txt", "w");
    for(int i = 0; i < steps; i++) {
        fprintf(fp5, "%.15f\t%.15f\t%.15f\n", x[i], y[i], tan(x[i]) - y[i]);
    }
    fclose(fp5);

    x2[0] = 0.0; // Reset initial x value
    y2[0] = 0.0; // Reset initial y value

    // Perform RK4 method with h = 0.05
    rk4(f, x2, y2, h2, steps2);

    FILE *fp6 = fopen("rk4_results_h_0.05.txt", "w");
    for(int i = 0; i < steps2; i++) {
        fprintf(fp6, "%.15f\t%.15f\t%.15f\n", x2[i], y2[i], tan(x2[i]) - y2[i]);
    }
    fclose(fp6);

    x3[0] = 0.0; // Reset initial x value
    y3[0] = 0.0; // Reset initial y value

    // Perform RK4 method with h = 0.01
    rk4(f, x3, y3, h3, steps3);

    FILE *fp7 = fopen("rk4_results_h_0.01.txt", "w");
    for(int i = 0; i < steps3; i++) {
        fprintf(fp7, "%.15f\t%.15f\t%.15f\n", x3[i], y3[i], tan(x3[i]) - y3[i]);
    }
    fclose(fp7);

    x4[0] = 0.0; // Reset initial x value
    y4[0] = 0.0; // Reset initial y value

    // Perform RK4 method with h = 0.005
    rk4(f, x4, y4, h4, steps4);

    FILE *fp8 = fopen("rk4_results_h_0.005.txt", "w");
    for(int i = 0; i < steps4; i++) {
        fprintf(fp8, "%.15f\t%.15f\t%.15f\n", x4[i], y4[i], tan(x4[i]) - y4[i]);
    }
    fclose(fp8);

    return 0;
}