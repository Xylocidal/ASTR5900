import numpy as np
import matplotlib.pyplot as plt

directory = '/Users/xylomolenda/Desktop/ASTR5900/Homework_01/'

data = np.genfromtxt(directory + 'HW01_data.txt', delimiter='\t', skip_header=1)
x = data[:, 0]
y = data[:, 1]

print("x:", x)
print("y:", y)