import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Constants
GM_sun = 4 * np.pi**2
GM_jupiter = GM_sun * 9.5e-4

# Acceleration functions
def accel_jupiter(r):
    r_mag = np.linalg.norm(r)
    return -GM_sun * r / r_mag**3

def accel_voyager(r_v, r_j):
    r_v_mag = np.linalg.norm(r_v)
    
    # Sun contribution
    a_sun = -GM_sun * r_v / r_v_mag**3
    
    # Jupiter contribution
    diff = r_v - r_j
    diff_mag = np.linalg.norm(diff)
    a_jup = -GM_jupiter * diff / diff_mag**3
    
    return a_sun + a_jup

# Leapfrog integrator (for both bodies)
def simulate(dt, steps, theta=np.pi / 2 + 0.179):
    # Jupiter initial conditions (circular orbit at 5.2 AU)
    r_j = np.array([5.2 * np.cos(theta), 5.2 * np.sin(theta)])
    v_j = np.array([np.sqrt(GM_sun / 5.2) * -np.sin(theta), np.sqrt(GM_sun / 5.2) * np.cos(theta)])
    
    # Voyager initial conditions (Earth position)
    r_v = np.array([1.0, 0.0])
    
    # Escape velocity from the Sun at 1 AU is sqrt(2*GM_sun/1.0)
    v_v = np.array([0.0, 0.95 * np.sqrt(2.0 * GM_sun / 1.0)])
    
    # Storage
    rj_hist = np.zeros((steps, 2))
    rv_hist = np.zeros((steps, 2))
    vv_hist = np.zeros((steps, 2))
    
    rj_hist[0] = r_j
    rv_hist[0] = r_v
    vv_hist[0] = v_v
    
    # Initial accelerations
    a_j = accel_jupiter(r_j)
    a_v = accel_voyager(r_v, r_j)
    
    # Half-step velocities
    vj_half = v_j + 0.5 * dt * a_j
    vv_half = v_v + 0.5 * dt * a_v
    
    for i in range(1, steps):
        r_j = r_j + dt * vj_half
        r_v = r_v + dt * vv_half
        
        a_j = accel_jupiter(r_j)
        a_v = accel_voyager(r_v, r_j)
        
        vj_half = vj_half + dt * a_j
        vv_half = vv_half + dt * a_v
        
        rj_hist[i] = r_j
        rv_hist[i] = r_v
        vv_hist[i] = vv_half - 0.5 * dt * a_v
    
    return rj_hist, rv_hist, vv_hist

# Simulation parameters
dt = 0.02
tmax = 6.0   # years
steps = int(tmax / dt)

r_j, r_v, v_v = simulate(dt, steps)

# Convert speed to km/s
AU_per_year_to_kms = 4.74047
speed = np.linalg.norm(v_v, axis=1) * AU_per_year_to_kms

# Animation
fig, ax = plt.subplots(figsize=(8, 8))

def update(i):
    ax.clear()
    
    ax.plot(r_j[:i, 0], r_j[:i, 1], 'orange', label='Jupiter')
    ax.plot(r_v[:i, 0], r_v[:i, 1], 'blue', label='Voyager 2')
    
    ax.scatter(r_j[i, 0], r_j[i, 1], color='orange')
    ax.scatter(r_v[i, 0], r_v[i, 1], color='blue')
    
    ax.scatter(0, 0, color='yellow', s=120, label='Sun')
    
    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    
    # Time in months
    time_months = i * dt * 12
    
    ax.set_title(f"Time = {time_months:.1f} months\nSpeed = {speed[i]:.2f} km/s")
    
    ax.legend()
    ax.grid()

ani = FuncAnimation(fig, update, frames=steps, interval=20)
ani.save("voyager.gif", writer=PillowWriter(fps=30))

# Plot speed and distance from Jupiter
distance_from_jupiter = np.linalg.norm(r_v - r_j, axis=1)
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(np.arange(steps) * dt * 12, speed, 'b-')
plt.xlabel('Time (months)')
plt.ylabel('Speed (km/s)')
plt.ylim(0, 40)
plt.title('Voyager 2 Speed')
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(np.arange(steps) * dt * 12, distance_from_jupiter, 'r-')
plt.xlabel('Time (months)')
plt.ylabel('Distance from Jupiter (AU)')
plt.title('Distance from Jupiter')
plt.grid()
plt.tight_layout()
plt.savefig("voyager_speed_distance.png")
plt.show()