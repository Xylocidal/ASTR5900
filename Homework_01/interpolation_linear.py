import numpy as np
import matplotlib.pyplot as plt

# Set font sizes for plots
smallfont = 12
mediumfont = 14
largefont = 16

# Load the data
directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'

datatoggle = 2 #change to 1 for HW01_data.txt, change to 2 for HW01_data2.txt

if datatoggle == 1:
    datafile = 'HW01_data.txt'
elif datatoggle == 2:
    datafile = 'HW01_data2.txt'

data = np.genfromtxt(directory + datafile, delimiter='\t', skip_header=1)
x = data[:, 0]
y = data[:, 1]

# Produce linear interpolation of data with given resolution
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
if datatoggle == 1:
    resolution = 81
elif datatoggle == 2:
    resolution = 91
x_linear, y_linear = linear_interpolation(x, y, resolution)

# Save the relative error to a text file
if datatoggle == 2:
    relative_error = np.abs((y_linear - (np.sin(0.5 * np.pi * x_linear) + 0.5 * x_linear)) / (np.sin(0.5 * np.pi * x_linear) + 0.5 * x_linear))
    error_data = np.column_stack((x_linear, relative_error))
    np.savetxt(directory + 'relative_error_linear' + str(datatoggle) + '.txt', error_data, delimiter='\t')

# Plot the original data and the linear interpolation
plt.figure(figsize=(5, 3))
plt.scatter(x, y, label='Original Data')
plt.scatter(x_linear, y_linear, s = 5, label='Linear Interpolation')

plt.title('Linear Interpolation of Data', fontsize = largefont)
plt.xlabel('x', fontsize = mediumfont)
plt.ylabel('y', fontsize = mediumfont)
plt.legend(fontsize = smallfont)
plt.grid()
plt.tight_layout()

plt.savefig(directory + 'interpolation_linear' + str(datatoggle) + '.png')