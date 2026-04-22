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

# estimate pi using Bayesian inference

# theta grid
theta = np.linspace(1e-6, 1-1e-6, 1000)

# prior (Beta(alpha,beta))
alpha, beta = 1, 1
prior = theta**(alpha-1) * (1-theta)**(beta-1)

# likelihood
H = np.sum(inside_circle)  # number of points inside the circle
T = len(inside_circle) - H  # number of points outside the circle
likelihood = theta**H * (1-theta)**T

# posterior
posterior_unnorm = prior * likelihood
posterior = posterior_unnorm / np.trapezoid(posterior_unnorm, theta)

# find the expected value of theta
theta_estimate = np.trapezoid(theta * posterior, theta)
pi_estimate = 4 * theta_estimate

print(f"Estimated value of pi: {pi_estimate}")