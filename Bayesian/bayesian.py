import numpy as np
import matplotlib.pyplot as plt

# load data
data = np.loadtxt("/Users/xylomolenda/Desktop/ASTR5900/Bayesian/HW06_data.txt")
H = np.sum(data)
T = len(data) - H

# define theta grid
theta = np.linspace(1e-6, 1-1e-6, 1000)

# define prior (arbitrary Beta distribution)
alpha, beta = 1, 1
prior = theta**(alpha-1) * (1-theta)**(beta-1)

# define likelihood function
likelihood = theta**H * (1-theta)**T

# calculate unnormalized posterior and normalize (Bayes' theorem)
posterior_unnorm = prior * likelihood
posterior = posterior_unnorm / np.trapezoid(posterior_unnorm, theta)

# plot the posterior distribution
plt.plot(theta, posterior)
plt.xlabel(r'$\theta$')
plt.ylabel('Posterior PDF')
plt.title('Posterior Distribution for Coin Bias')
plt.savefig("posterior_distribution.png")
plt.show()