#include <stdio.h>
#include <math.h>

double f(double x, double y);
double probability_density(double m, double T, double v);
void rk4(double (*f)(double, double), double *x, double *y, double h, int steps);
void euler(double (*f)(double, double), double *x, double *y, double h, int steps);

typedef struct {
    double m;
    double T;
} PDFParams;

double pdf_rhs(double v, double I, void *params);

void rk4_ctx(double (*f)(double, double, void *),
             void *params,
             double *x, double *y,
             double h, int steps);

double f(double x, double y) {
    return y * y + 1.0;
}

double probability_density(double m, double T, double v) {
    const double k_B = 1.380649e-23; // Boltzmann constant in J/K
    double exponent = -0.5 * m * v * v / (k_B * T);
    return pow(m / (2.0 * M_PI * k_B * T), 1.5) * 4 * M_PI * v * v * exp(exponent);
}

double pdf_rhs(double v, double I, void *params) {
    PDFParams *p = (PDFParams *)params;
    return probability_density(p->m, p->T, v);
}

void rk4_ctx(double (*f)(double, double, void *),
             void *params,
             double *x, double *y,
             double h, int steps) {
    double k1, k2, k3, k4;

    for (int i = 0; i < steps - 1; i++) {
        k1 = f(x[i],           y[i],           params);
        k2 = f(x[i] + h / 2.0, y[i] + h*k1/2.0, params);
        k3 = f(x[i] + h / 2.0, y[i] + h*k2/2.0, params);
        k4 = f(x[i] + h,       y[i] + h*k3,     params);

        y[i + 1] = y[i] + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        x[i + 1] = x[i] + h;
    }
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

    // Print probability density to a file over the interval [0, 1000] m/s for a Hydrogen atom at 10000K
    double m = 1.6735575e-27; // Mass of a Hydrogen atom in kg
    double T = 10000.0; // Temperature in K
    FILE *fp9 = fopen("probability_density.txt", "w");
    for (double v = 0.0; v <= 50000.0; v += 50.0) {
        fprintf(fp9, "%.2f\t%.15e\n", v, probability_density(m, T, v));
    }
    fclose(fp9);

    // Integrate the probability density over the speed range using RK4:
    // I'(v) = p(v), I(0)=0. Then I(v_max) should be ~ 1.
    int steps5 = 1001;          // include endpoint nicely
    double v_grid[steps5];
    double integral[steps5];
    double hstar = 50.0;        // m/s
    double v_max = hstar * (steps5 - 1);

    v_grid[0] = 0.0;
    integral[0] = 0.0;

    PDFParams params = { .m = m, .T = T };

    rk4_ctx(pdf_rhs, &params, v_grid, integral, hstar, steps5);

    printf("Integral from 0 to %.1f m/s = %.15f\n", v_max, integral[steps5 - 1]);

    return 0;
}