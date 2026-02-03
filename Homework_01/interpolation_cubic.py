import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

# Set font sizes for plots
smallfont = 12
mediumfont = 14
largefont = 16

# Loads the data
directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'

data = np.genfromtxt(directory + 'HW01_data.txt', delimiter='\t', skip_header=1)
x = data[:, 0]
y = data[:, 1]

# Perform cubic spline interpolation on the given data using a scipy library function
cubic_spline = interp.CubicSpline(x, y)

# Generate new x values for plotting the cubic spline
x_spline = np.linspace(x.min(), x.max(), 100)
y_spline = cubic_spline(x_spline)

# Plot the original data and the cubic spline interpolation
plt.figure(figsize=(10, 6))
plt.scatter(x, y, label='Original Data')
plt.plot(x_spline, y_spline, label='Cubic Spline Interpolation')

plt.title('Cubic Spline Interpolation of Data', fontsize = largefont)
plt.xlabel('x', fontsize = mediumfont)
plt.ylabel('y', fontsize = mediumfont)
plt.legend(fontsize = smallfont)
plt.grid()

plt.savefig(directory + 'interpolation_cubic.png')