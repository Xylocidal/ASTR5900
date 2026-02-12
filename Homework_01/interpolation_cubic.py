import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

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

# Perform cubic spline interpolation on the given data using a scipy library function
cubic_spline = interp.CubicSpline(x, y)

# Generate new x values for plotting the cubic spline
if datatoggle == 1:
    resolution = 81
elif datatoggle == 2:
    resolution = 101
x_spline = np.linspace(x.min(), x.max(), resolution)
y_spline = cubic_spline(x_spline)

# Save the relative error to a text file
if datatoggle == 2:
    relative_error = np.zeros(resolution)
    for i in range(resolution):
        if i == 0:
            relative_error[i] = 0 # Set relative error to 0 at the first point to avoid division by zero
        else:
            relative_error[i] = (y_spline[i] - (np.sin(0.5 * np.pi * x_spline[i]) + 0.5 * x_spline[i])) / (np.sin(0.5 * np.pi * x_spline[i]) + 0.5 * x_spline[i])
    error_data = np.column_stack((x_spline, relative_error))
    np.savetxt(directory + 'relative_error_cubic' + str(datatoggle) + '.txt', error_data, delimiter='\t')

# Plot the original data and the cubic spline interpolation
plt.figure(figsize=(5, 3))
plt.scatter(x, y, label='Original Data')
plt.plot(x_spline, y_spline, label='Cubic Spline Interpolation')

plt.title('Cubic Spline Interpolation of Data', fontsize = largefont)
plt.xlabel('x', fontsize = mediumfont)
plt.ylabel('y', fontsize = mediumfont)
plt.legend(fontsize = smallfont)
plt.grid()
plt.tight_layout()

plt.savefig(directory + 'interpolation_cubic' + str(datatoggle) + '.png')