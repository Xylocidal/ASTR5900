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