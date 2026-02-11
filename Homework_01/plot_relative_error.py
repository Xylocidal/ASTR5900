import numpy as np
import matplotlib.pyplot as plt

# Set font sizes for plots
smallfont = 12
mediumfont = 14
largefont = 16

# Load the data
directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'

linearerror = np.genfromtxt(directory + "relative_error_linear2.txt", delimiter='\t')
cubicerror = np.genfromtxt(directory + "relative_error_cubic2.txt", delimiter='\t')

# Plot the relative error for linear and cubic interpolation
plt.figure(figsize=(5, 3))
plt.scatter(cubicerror[:, 0], cubicerror[:, 1], label='Cubic Spline', s = 5)
plt.scatter(linearerror[:, 0], linearerror[:, 1], label='Linear', s = 5)


plt.title('Relative Error of Interpolation Methods', fontsize = largefont)
plt.xlabel('x', fontsize = mediumfont)
plt.ylabel('Relative Error', fontsize = mediumfont)
plt.legend(fontsize = smallfont)
plt.grid()
plt.tight_layout()
plt.savefig(directory + 'relative_error_comparison.png')
