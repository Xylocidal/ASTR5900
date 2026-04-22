import numpy as np
import matplotlib.pyplot as plt

# load data
data = np.loadtxt("/Users/xylomolenda/Desktop/ASTR5900/Bayesian/HW06_data.txt")

# theta grid
theta = np.linspace(1e-6, 1-1e-6, 1000)

# prior (Beta(alpha,beta))
alpha, beta = 5, 10
prior = theta**(alpha-1) * (1-theta)**(beta-1)

# values of N to compare
Ns = [5, 50, 500]

# linestyles for plotting
linestyles = ['-', '--', ':']

plt.figure()

for N in Ns:
    H = np.sum(data[:N])
    T = N - H

    likelihood = theta**H * (1-theta)**T
    posterior_unnorm = prior * likelihood
    posterior = posterior_unnorm / np.trapezoid(posterior_unnorm, theta)

    plt.plot(theta, posterior, label=f"N={N}", linestyle=linestyles[Ns.index(N)])

# plot formatting
plt.xlabel(r'$\theta$', fontsize=14)
plt.ylabel('Posterior PDF', fontsize=14)
plt.title('Posterior Distribution for Coin Bias with Beta(5,10) Prior', fontsize=14)
plt.legend(fontsize=14)

plt.savefig("posterior_comparison_beta5_10.png")
plt.show()