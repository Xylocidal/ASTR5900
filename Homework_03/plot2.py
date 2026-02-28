import numpy as np
import matplotlib.pyplot as plt

# Font sizes
smallfont = 12
mediumfont = 14
largefont = 16

# Load data from files
directory = "/Users/xylomolenda/Desktop/ASTR5900/Homework_03/"
euler_data_h_0_1   = np.loadtxt(directory + "euler_results_h_0.1.txt")
euler_data_h_0_05  = np.loadtxt(directory + "euler_results_h_0.05.txt")
euler_data_h_0_01  = np.loadtxt(directory + "euler_results_h_0.01.txt")
euler_data_h_0_005 = np.loadtxt(directory + "euler_results_h_0.005.txt")
rk4_data_h_0_1     = np.loadtxt(directory + "rk4_results_h_0.1.txt")
rk4_data_h_0_05    = np.loadtxt(directory + "rk4_results_h_0.05.txt")
rk4_data_h_0_01    = np.loadtxt(directory + "rk4_results_h_0.01.txt")
rk4_data_h_0_005   = np.loadtxt(directory + "rk4_results_h_0.005.txt")

# Extract x, y, and error values
x_euler_h_0_1 = euler_data_h_0_1[:, 0]
y_euler_h_0_1 = euler_data_h_0_1[:, 1]
error_euler_h_0_1 = euler_data_h_0_1[:, 2]

x_euler_h_0_05 = euler_data_h_0_05[:, 0]
y_euler_h_0_05 = euler_data_h_0_05[:, 1]
error_euler_h_0_05 = euler_data_h_0_05[:, 2]

x_euler_h_0_01 = euler_data_h_0_01[:, 0]
y_euler_h_0_01 = euler_data_h_0_01[:, 1]
error_euler_h_0_01 = euler_data_h_0_01[:, 2]

x_euler_h_0_005 = euler_data_h_0_005[:, 0]
y_euler_h_0_005 = euler_data_h_0_005[:, 1]
error_euler_h_0_005 = euler_data_h_0_005[:, 2]

x_rk4_h_0_1 = rk4_data_h_0_1[:, 0]
y_rk4_h_0_1 = rk4_data_h_0_1[:, 1]
error_rk4_h_0_1 = rk4_data_h_0_1[:, 2]

x_rk4_h_0_05 = rk4_data_h_0_05[:, 0]
y_rk4_h_0_05 = rk4_data_h_0_05[:, 1]
error_rk4_h_0_05 = rk4_data_h_0_05[:, 2]

x_rk4_h_0_01 = rk4_data_h_0_01[:, 0]
y_rk4_h_0_01 = rk4_data_h_0_01[:, 1]
error_rk4_h_0_01 = rk4_data_h_0_01[:, 2]

x_rk4_h_0_005 = rk4_data_h_0_005[:, 0]
y_rk4_h_0_005 = rk4_data_h_0_005[:, 1]
error_rk4_h_0_005 = rk4_data_h_0_005[:, 2]


def fractional_solution_error_vs_reference(x_coarse, y_coarse, x_ref, y_ref, eps=1e-15):
    # Compute |y_coarse(x) - y_ref(x)| / |y_ref(x)|,
    # where y_ref(x) is interpolated onto the coarse x-grid.
    # eps prevents division by zero when y_ref ~ 0.
    y_ref_on_coarse = np.interp(x_coarse, x_ref, y_ref)
    denom = np.maximum(np.abs(y_ref_on_coarse), eps)
    return np.abs(y_coarse - y_ref_on_coarse) / denom


# Fractional solution error relative to h=0.005 reference, evaluated on each coarse grid
fractional_error_euler_h_0_1  = fractional_solution_error_vs_reference(
    x_euler_h_0_1,  y_euler_h_0_1,  x_euler_h_0_005, y_euler_h_0_005
)
fractional_error_euler_h_0_05 = fractional_solution_error_vs_reference(
    x_euler_h_0_05, y_euler_h_0_05, x_euler_h_0_005, y_euler_h_0_005
)
fractional_error_euler_h_0_01 = fractional_solution_error_vs_reference(
    x_euler_h_0_01, y_euler_h_0_01, x_euler_h_0_005, y_euler_h_0_005
)

fractional_error_rk4_h_0_1  = fractional_solution_error_vs_reference(
    x_rk4_h_0_1,  y_rk4_h_0_1,  x_rk4_h_0_005, y_rk4_h_0_005
)
fractional_error_rk4_h_0_05 = fractional_solution_error_vs_reference(
    x_rk4_h_0_05, y_rk4_h_0_05, x_rk4_h_0_005, y_rk4_h_0_005
)
fractional_error_rk4_h_0_01 = fractional_solution_error_vs_reference(
    x_rk4_h_0_01, y_rk4_h_0_01, x_rk4_h_0_005, y_rk4_h_0_005
)


# Create plots

# Euler: |y_h - y_ref| / |y_ref|
plt.figure(figsize=(12, 6))
plt.plot(x_euler_h_0_1,  fractional_error_euler_h_0_1,  label=r"$|y_{0.1}-y_{0.005}|/|y_{0.005}|$",  marker='o')
plt.plot(x_euler_h_0_05, fractional_error_euler_h_0_05, label=r"$|y_{0.05}-y_{0.005}|/|y_{0.005}|$", marker='o')
plt.plot(x_euler_h_0_01, fractional_error_euler_h_0_01, label=r"$|y_{0.01}-y_{0.005}|/|y_{0.005}|$", marker='o')

# Plot reference error relative to exact solution
plt.plot(x_euler_h_0_005, np.abs(error_euler_h_0_005), label=r"$|y_{0.005}-y_{\rm exact}|/|y_{\rm exact}|$", marker='o')

plt.yscale('log')
plt.xscale('log')
plt.ylim(1e-9, 1e3)
plt.title("Fractional Error (Euler)", fontsize=largefont)
plt.xlabel("x", fontsize=mediumfont)
plt.ylabel("Error", fontsize=mediumfont)
plt.legend(fontsize=smallfont)
plt.grid()
plt.savefig(directory + "fractional_error_euler.png", dpi=300, bbox_inches="tight")


# RK4: |y_h - y_ref| / |y_ref|
plt.figure(figsize=(12, 6))
plt.plot(x_rk4_h_0_1,  fractional_error_rk4_h_0_1,  label=r"$|y_{0.1}-y_{0.005}|/|y_{0.005}|$",  marker='o')
plt.plot(x_rk4_h_0_05, fractional_error_rk4_h_0_05, label=r"$|y_{0.05}-y_{0.005}|/|y_{0.005}|$", marker='o')
plt.plot(x_rk4_h_0_01, fractional_error_rk4_h_0_01, label=r"$|y_{0.01}-y_{0.005}|/|y_{0.005}|$", marker='o')

# Plot reference error relative to exact solution
plt.plot(x_rk4_h_0_005, np.abs(error_rk4_h_0_005), label=r"$|y_{0.005}-y_{\rm exact}|/|y_{\rm exact}|$", marker='o')

plt.yscale('log')
plt.xscale('log')
plt.ylim(1e-16, 1e3)
plt.title("Fractional Error (RK4)", fontsize=largefont)
plt.xlabel("x", fontsize=mediumfont)
plt.ylabel("Error", fontsize=mediumfont)
plt.legend(fontsize=smallfont)
plt.grid()
plt.savefig(directory + "fractional_error_rk4.png", dpi=300, bbox_inches="tight")