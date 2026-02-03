import numpy as np
import matplotlib.pyplot as plt

# Loads the data
directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'

data = np.genfromtxt(directory + 'HW01_data.txt', delimiter='\t', skip_header=1)
x = data[:, 0]
y = data[:, 1]

# Produces linear interpolation of data with given resolution
def linear_interpolation(x, y, resolution):
    n = len(x)
    b = y[:-1] #b_i coefficients
    a = np.zeros(n - 1)

    #calculates a_i coefficients
    for i in range(n - 1):
        a[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])

    #generates new x values in original interval based on resolution
    x_new = np.linspace(x[0], x[-1], resolution)
    y_new = np.zeros(resolution)

    #calculates y_new values based on linear interpolation
    for i in range(n - 1):
        for j in range(resolution):
            if x_new[j] >= x[i] and x_new[j] <= x[i + 1]:
                y_new[j] = a[i] * (x_new[j] - x[i]) + b[i]
    return x_new, y_new

# Perform linear interpolation on the given data
resolution = 100
x_linear, y_linear = linear_interpolation(x, y, resolution)
