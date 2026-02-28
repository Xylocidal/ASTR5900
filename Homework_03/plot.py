import numpy as np
import matplotlib.pyplot as plt

# Font sizes
smallfont = 16
mediumfont = 18
largefont = 20

# Calculate tan(x) for plotting the exact solution
x_values = np.linspace(0, 2, 1000)
tan_values = np.tan(x_values)

# Load data from files
directory = "/Users/xylomolenda/Desktop/ASTR5900/Homework_03/"
euler_data_h_0_1 = np.loadtxt(directory + "euler_results_h_0.1.txt")
rk4_data_h_0_1 = np.loadtxt(directory + "rk4_results_h_0.1.txt")
euler_data_h_0_05 = np.loadtxt(directory + "euler_results_h_0.05.txt")
rk4_data_h_0_05 = np.loadtxt(directory + "rk4_results_h_0.05.txt")

# Extract x, y, and error values
x_euler_h_0_1 = euler_data_h_0_1[:, 0]
y_euler_h_0_1 = euler_data_h_0_1[:, 1]
error_euler_h_0_1 = euler_data_h_0_1[:, 2] 

x_rk4_h_0_1 = rk4_data_h_0_1[:, 0]
y_rk4_h_0_1 = rk4_data_h_0_1[:, 1]
error_rk4_h_0_1 = rk4_data_h_0_1[:, 2]

x_euler_h_0_05 = euler_data_h_0_05[:, 0]
y_euler_h_0_05 = euler_data_h_0_05[:, 1]
error_euler_h_0_05 = euler_data_h_0_05[:, 2]

x_rk4_h_0_05 = rk4_data_h_0_05[:, 0]
y_rk4_h_0_05 = rk4_data_h_0_05[:, 1]
error_rk4_h_0_05 = rk4_data_h_0_05[:, 2]

# Create plots

# Plot h = 0.1 data
plt.figure(figsize=(12, 6))
plt.plot(x_values, tan_values, label="tan(x)", linestyle='-')
plt.plot(x_euler_h_0_1, y_euler_h_0_1, label="Euler", marker='o')
plt.plot(x_rk4_h_0_1, y_rk4_h_0_1, label="RK4", marker='D')
plt.ylim(-10, 200)  # Limit y-axis to avoid extreme values
plt.title("Numerical Integration (h = 0.1)", fontsize=largefont)
plt.xlabel("x", fontsize=mediumfont)
plt.ylabel("y", fontsize=mediumfont)
plt.xticks(fontsize=smallfont)
plt.yticks(fontsize=smallfont)
plt.legend(fontsize=smallfont)
plt.grid()
plt.savefig(directory + "numerical_integration_h_0_1.png")

# Plot h = 0.05 data
plt.figure(figsize=(12, 6))
plt.plot(x_values, tan_values, label="tan(x)", linestyle='-')
plt.plot(x_euler_h_0_05, y_euler_h_0_05, label="Euler", marker='o')
plt.plot(x_rk4_h_0_05, y_rk4_h_0_05, label="RK4", marker='D')
plt.ylim(-10, 200)  # Limit y-axis to avoid extreme values
plt.title("Numerical Integration (h = 0.05)", fontsize=largefont)
plt.xlabel("x", fontsize=mediumfont)
plt.ylabel("y", fontsize=mediumfont)
plt.xticks(fontsize=smallfont)
plt.yticks(fontsize=smallfont)
plt.legend(fontsize=smallfont)
plt.grid()
plt.savefig(directory + "numerical_integration_h_0_05.png")