import numpy as np
import matplotlib.pyplot as plt

# Diffusion constant
alpha = 0.01

# Initial condition (periodic on [0,1])
def u0(x):
    return np.exp(-50*(x-0.5)**2)

# Discrete Fourier Transform
def dft(x):
    N = len(x)
    X = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            X[k] += x[n] * np.exp(-2j * np.pi * k * n / N)
    return X

# Inverse Discrete Fourier Transform
def idft(X):
    N = len(X)
    x = np.zeros(N, dtype=complex)
    for n in range(N):
        for k in range(N):
            x[n] += X[k] * np.exp(2j * np.pi * k * n / N)
    return x / N

# RK4 method for solving du/dt = -alpha * (2*pi*k)^2 * u in Fourier space
def rk4_step(u_hat, dt, k_vals):

    def rhs(u):
        return -alpha * (2*np.pi*k_vals)**2 * u

    k1 = rhs(u_hat)
    k2 = rhs(u_hat + 0.5*dt*k1)
    k3 = rhs(u_hat + 0.5*dt*k2)
    k4 = rhs(u_hat + dt*k3)

    return u_hat + dt*(k1 + 2*k2 + 2*k3 + k4)/6


# Spatial grid
N = 128
x = np.linspace(0, 1, N, endpoint=False)

# Time parameters
dt = 0.001
tmax = 10.0
steps = int(tmax/dt)

# Initial condition
u = u0(x)

# Fourier transform initial condition
u_hat = dft(u)

# Fourier mode indices
k_vals = np.arange(N)
k_vals[k_vals > N//2] -= N # shift to negative frequencies

# Storage array
u_xt = np.zeros((steps, N))

# Time evolution
for i in range(steps):

    u = np.real(idft(u_hat)) # transform back to real space at this step
    u_xt[i] = u

    u_hat = rk4_step(u_hat, dt, k_vals) # evolve to next step in Fourier space

# Time axis
t = np.linspace(0, tmax, steps)

# Create meshgrid
X, T = np.meshgrid(x, t)

# Plot using pcolormesh
plt.figure(figsize=(8,6))

pcm = plt.pcolormesh(
    X, T, u_xt,
    shading='auto'
)

plt.xlabel("x")
plt.ylabel("t")
plt.title("Heat Equation Solution")

plt.colorbar(pcm, label="u(x,t)")

plt.savefig("heat_equation_density.png")
plt.show()

print("heat values at t=10:", u_xt[-1])