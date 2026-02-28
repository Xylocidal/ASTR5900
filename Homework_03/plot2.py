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
euler_data_h_0_001 = np.loadtxt(directory + "euler_results_h_0.001.txt")

rk4_data_h_0_1     = np.loadtxt(directory + "rk4_results_h_0.1.txt")
rk4_data_h_0_05    = np.loadtxt(directory + "rk4_results_h_0.05.txt")
rk4_data_h_0_01    = np.loadtxt(directory + "rk4_results_h_0.01.txt")
rk4_data_h_0_005   = np.loadtxt(directory + "rk4_results_h_0.005.txt")
rk4_data_h_0_001   = np.loadtxt(directory + "rk4_results_h_0.001.txt")

# Extract x and y
x_euler_0_1,   y_euler_0_1   = euler_data_h_0_1[:,0],   euler_data_h_0_1[:,1]
x_euler_0_05,  y_euler_0_05  = euler_data_h_0_05[:,0],  euler_data_h_0_05[:,1]
x_euler_0_01,  y_euler_0_01  = euler_data_h_0_01[:,0],  euler_data_h_0_01[:,1]
x_euler_0_005, y_euler_0_005 = euler_data_h_0_005[:,0], euler_data_h_0_005[:,1]
x_euler_0_001, y_euler_0_001 = euler_data_h_0_001[:,0], euler_data_h_0_001[:,1]

x_rk4_0_1,   y_rk4_0_1   = rk4_data_h_0_1[:,0],   rk4_data_h_0_1[:,1]
x_rk4_0_05,  y_rk4_0_05  = rk4_data_h_0_05[:,0],  rk4_data_h_0_05[:,1]
x_rk4_0_01,  y_rk4_0_01  = rk4_data_h_0_01[:,0],  rk4_data_h_0_01[:,1]
x_rk4_0_005, y_rk4_0_005 = rk4_data_h_0_005[:,0], rk4_data_h_0_005[:,1]
x_rk4_0_001, y_rk4_0_001 = rk4_data_h_0_001[:,0], rk4_data_h_0_001[:,1]

#linear interpolation function to get y at any x from the grid of (x,y) points
def y_at_x(x_query, x_grid, y_grid):
    return np.interp(x_query, x_grid, y_grid)

# function to compute fractional difference at a specific x_query
def frac_diff_at_x(x_query, x_coarse, y_coarse, x_ref, y_ref, eps=1e-15):
    yc = y_at_x(x_query, x_coarse, y_coarse)
    yr = y_at_x(x_query, x_ref, y_ref)
    return np.abs(yc - yr) / max(np.abs(yr), eps)


x_f = 1.0   # x used for convergence study
eps = 1e-15 # small number to avoid division by zero in fractional difference

# Step sizes and corresponding datasets (coarse -> reference is h=0.005)
hs = np.array([0.1, 0.05, 0.01, 0.005])
# Reference for each method
x_ref_e, y_ref_e = x_euler_0_001, y_euler_0_001
x_ref_r, y_ref_r = x_rk4_0_001,   y_rk4_0_001

# Compute fractional differences at x_f for each h
err_euler = np.array([
    frac_diff_at_x(x_f, x_euler_0_1,  y_euler_0_1,  x_ref_e, y_ref_e, eps=eps),
    frac_diff_at_x(x_f, x_euler_0_05, y_euler_0_05, x_ref_e, y_ref_e, eps=eps),
    frac_diff_at_x(x_f, x_euler_0_01, y_euler_0_01, x_ref_e, y_ref_e, eps=eps),
    frac_diff_at_x(x_f, x_euler_0_005, y_euler_0_005, x_ref_e, y_ref_e, eps=eps),
])

err_rk4 = np.array([
    frac_diff_at_x(x_f, x_rk4_0_1,  y_rk4_0_1,  x_ref_r, y_ref_r, eps=eps),
    frac_diff_at_x(x_f, x_rk4_0_05, y_rk4_0_05, x_ref_r, y_ref_r, eps=eps),
    frac_diff_at_x(x_f, x_rk4_0_01, y_rk4_0_01, x_ref_r, y_ref_r, eps=eps),
    frac_diff_at_x(x_f, x_rk4_0_005, y_rk4_0_005, x_ref_r, y_ref_r, eps=eps),
])

# Estimate convergence rates (slope in log-log space)
rate_euler = np.log(err_euler[:-1] / err_euler[1:]) / np.log(hs[:-1] / hs[1:])
rate_rk4 = np.log(err_rk4[:-1] / err_rk4[1:]) / np.log(hs[:-1] / hs[1:])

print("Estimated convergence rates:")
for i in range(len(hs)-1):
    print(f"Euler: h={hs[i]:.3f} to h={hs[i+1]:.3f} -> rate ~ {rate_euler[i]:.2f}")
    print(f"RK4:   h={hs[i]:.3f} to h={hs[i+1]:.3f} -> rate ~ {rate_rk4[i]:.2f}")

# Plotting the convergence study on a log-log scale
plt.figure(figsize=(8,6))
plt.loglog(hs, err_euler, marker='o', label="Euler", color = "tab:orange")
plt.loglog(hs, err_rk4,   marker='D', label="RK4", color = "tab:green")

plt.title(r"Convergence study: fractional difference vs step size $h$" + f"\n(evaluated at x = {x_f})",
          fontsize=largefont)
plt.xlabel(r"Step size $h$", fontsize=mediumfont)
plt.ylabel(r"Fractional difference $|y_h - y_{\rm ref}|/|y_{\rm ref}|$", fontsize=mediumfont)
plt.grid(True, which="both")
plt.legend(fontsize=smallfont)

plt.savefig(directory + "convergence_loglog_euler_vs_rk4.png", dpi=300, bbox_inches="tight")
plt.show()