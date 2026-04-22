import numpy as np
import matplotlib.pyplot as plt

np.random.seed(4148379212)

Ns = [10, 100, 1000]
pi_estimates = {}
x_full = np.random.rand(1000)
y_full = np.random.rand(1000)

# theta grid for Bayesian step
theta = np.linspace(1e-6, 1-1e-6, 1000)

alpha, beta = 1, 1
prior = theta**(alpha-1) * (1-theta)**(beta-1)

for N in Ns:
    x = x_full[:N]
    y = y_full[:N]

    inside = x**2 + y**2 <= 1

    # only plot for N=100
    if N == 100:
        plt.figure(figsize=(6, 6))
        plt.scatter(x[inside], y[inside], color='blue', label='Inside Circle')
        plt.scatter(x[~inside], y[~inside], color='red', label='Outside Circle')
        plt.xlabel('x', fontsize=14)
        plt.ylabel('y', fontsize=14)
        plt.grid()
        plt.savefig("random_points_circle.png")
        plt.show()

    # Bayesian estimate
    H = np.sum(inside)
    T = N - H

    likelihood = theta**H * (1-theta)**T
    posterior_unnorm = prior * likelihood
    posterior = posterior_unnorm / np.trapezoid(posterior_unnorm, theta)

    theta_est = np.trapezoid(theta * posterior, theta)
    pi_estimates[N] = 4 * theta_est

# print results
for N in Ns:
    print(f"N = {N}, estimated pi = {pi_estimates[N]:.6f}")