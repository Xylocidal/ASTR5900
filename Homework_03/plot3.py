import numpy as np
import matplotlib.pyplot as plt

# Font sizes
smallfont = 12
mediumfont = 14
largefont = 16

# Load data from files
directory = "/Users/xylomolenda/Desktop/ASTR5900/Homework_03/"
probability_density_data = np.loadtxt(directory + "probability_density.txt")

# Extract speed and probability density values
v = probability_density_data[:, 0]
probability_density = probability_density_data[:, 1]

# Create plots

# Plot probability density
plt.figure(figsize=(12, 6))
plt.plot(v, probability_density, label="Probability Density")
plt.title("Probability Density of Speeds for Hydrogen at 10000K", fontsize=largefont)
plt.xlabel("Speed (m/s)", fontsize=mediumfont)
plt.ylabel("Probability Density", fontsize=mediumfont)
plt.legend(fontsize=smallfont)
plt.grid()
plt.savefig(directory + "probability_density.png")