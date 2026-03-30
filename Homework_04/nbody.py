import numpy as np
import matplotlib.pyplot as plt

def euler_step(r, v, dt):
    """Perform one step of Euler's method."""
    r_mag = np.linalg.norm(r)
    a = -r / r_mag**3  # Gravitational acceleration
    r_new = r + v * dt
    v_new = v + a * dt
    return r_new, v_new

def leapfrog_step(r, v, dt):
    """Perform one step of the leapfrog method."""
    r_mag = np.linalg.norm(r)
    a = -r / r_mag**3  # Gravitational acceleration
    v_new = v + a * dt
    r_new = r + v_new * dt
    return r_new, v_new


# Initial conditions
r0 = np.array([1.0, 0.0])  # Initial position
v0 = np.array([0.0, 2.0 * np.pi])  # Initial velocity
dt = 0.01  # Time step
tmax = 3.0  # Total simulation time
steps = int(tmax / dt)

# Arrays to store the trajectory
r_euler = np.zeros((steps, 2))
v_euler = np.zeros((steps, 2))
r_leapfrog = np.zeros((steps, 2))
v_leapfrog = np.zeros((steps, 2))

# Set initial conditions
r_euler[0] = r0
v_euler[0] = v0
r_leapfrog[0] = r0
v_leapfrog[0] = v0

# Run the simulation
for i in range(1, steps):
    r_euler[i], v_euler[i] = euler_step(r_euler[i-1], v_euler[i-1], dt)
    r_leapfrog[i], v_leapfrog[i] = leapfrog_step(r_leapfrog[i-1], v_leapfrog[i-1], dt)

# Compute the specific total energy for both methods
def specific_energy(r, v):
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    return 0.5 * v_mag**2 - (4 * np.pi**2) / r_mag