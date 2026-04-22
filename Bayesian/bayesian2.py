import numpy as np
import matplotlib.pyplot as plt

# generate 100 random data points in first quadrant
np.random.seed(4148379210)
x = np.random.rand(100)
y = np.random.rand(100)

# calculate whether each point is inside the unit circle
inside_circle = x**2 + y**2 <= 1

# plot the data
plt.figure(figsize=(6, 6))
plt.scatter(x[inside_circle], y[inside_circle], color='blue', label='Inside Circle')
plt.scatter(x[~inside_circle], y[~inside_circle], color='red', label='Outside Circle')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.savefig("random_points_circle.png")
