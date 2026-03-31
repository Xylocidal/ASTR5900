import numpy as np
import matplotlib.pyplot as plt

def euler_step(r, v, dt):
    """Perform one step of Euler's method."""
    r_mag = np.linalg.norm(r)
    a = - (4 * np.pi**2) * r / r_mag**3  # Gravitational acceleration
    r_new = r + v * dt
    v_new = v + a * dt
    return r_new, v_new

def leapfrog(r, v, dt, steps):
    """Perform leapfrog integration."""
    r_out = np.zeros((steps, 2))
    v_out = np.zeros((steps, 2))
    
    r_out[0] = r
    v_out[0] = v
    
    # Initial acceleration
    r_mag = np.linalg.norm(r)
    a = - (4 * np.pi**2) * r / r_mag**3
    
    # Initialize half-step velocity
    v_half = v + 0.5 * dt * a

    for i in range(1, steps):
        r = r + dt * v_half
        r_mag = np.linalg.norm(r)
        a = - (4 * np.pi**2) * r / r_mag**3
        
        v_half = v_half + dt * a
        
        r_out[i] = r
        v_out[i] = v_half - 0.5 * dt * a

    return r_out, v_out


# Initial conditions
r0 = np.array([1.0, 0.0])  # Initial position
v0 = np.array([0.0, 0.8 * 2.0 * np.pi])  # Initial velocity
dt = 0.01  # Time step
tmax = 5.0  # Total simulation time
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

r_leapfrog, v_leapfrog = leapfrog(r_leapfrog[0], v_leapfrog[0], dt, steps)

# Plot the trajectories as a movie
from matplotlib.animation import FuncAnimation, PillowWriter

fig, ax = plt.subplots(figsize=(8, 8))

def update(i):
    ax.clear()
    
    ax.plot(r_euler[:i, 0], r_euler[:i, 1], 'r-', label='Euler')
    ax.plot(r_leapfrog[:i, 0], r_leapfrog[:i, 1], 'b-', label='Leapfrog')
    ax.scatter(r_euler[i, 0], r_euler[i, 1], color='red')
    ax.scatter(r_leapfrog[i, 0], r_leapfrog[i, 1], color='blue')
    
    ax.scatter(0, 0, color='yellow', s=100, label='Sun')
    
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Time = {i*dt:.2f} years')
    ax.legend()
    ax.grid()

ani = FuncAnimation(fig, update, frames=steps, interval=20)

# Save as GIF
ani.save("orbit2.gif", writer=PillowWriter(fps=30))

# Compute the specific total energy for both methods
def specific_energy(r, v):
    r_mag = np.linalg.norm(r, axis=1)
    v_mag = np.linalg.norm(v, axis=1)
    return 0.5 * v_mag**2 - (4 * np.pi**2) / r_mag

energy_euler = specific_energy(r_euler, v_euler)
energy_leapfrog = specific_energy(r_leapfrog, v_leapfrog)

# Plot the energy as a function of time
smallfont = 12
mediumfont = 14
largefont = 16

plt.figure(figsize=(10, 5))
plt.plot(np.arange(steps) * dt, energy_euler, 'r-', label='Euler Energy')
plt.plot(np.arange(steps) * dt, energy_leapfrog, 'b-', label='Leapfrog Energy')
plt.xlabel('Time (years)', fontsize=mediumfont)
plt.ylabel('Specific Total Energy (AU^2/year^2)', fontsize=mediumfont)
plt.title('Energy vs Time', fontsize=largefont)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("energy_speed2.png")
plt.show()