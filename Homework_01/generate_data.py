import numpy as np

def func(x):
    return np.sin(0.5 * np.pi * x) + 0.5 * x

# Generate x values from 0 to 10 with a step of 1
x = np.linspace(0, 10, 11)
y = func(x)

# Save the data to a text file
directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'
data = np.column_stack((x, y))
np.savetxt(directory + 'HW01_data2.txt', data, delimiter='\t')