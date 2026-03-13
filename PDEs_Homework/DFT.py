import numpy as np
import matplotlib.pyplot as plt
import time

# Input function f(x) Gaussian on [0, 1]
def f(x):
    return np.exp(-50*(x-0.5)**2)

# Analytical Fourier Transform of f(x)
def F(k):
    return np.sqrt(np.pi/50) * np.exp(- (np.pi**2 * k**2) / 50) * np.exp(-1j * np.pi * k)

# Discrete Fourier Transform
def dft(x):
    N = len(x)
    X = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            X[k] += 1/N * x[n] * np.exp(-1j * 2 * np.pi * k * n / N)
    return X

# Inverse Discrete Fourier Transform
def idft(X):
    N = len(X)
    x = np.zeros(N, dtype=complex)
    for n in range(N):
        for k in range(N):
            x[n] +=  N * X[k] * np.exp(1j * 2 * np.pi * k * n / N)
    return x / N

# Sample the function at three different values of N
N1 = 16
x_samples1 = np.linspace(0, 1, N1, endpoint=False)
f_samples1 = f(x_samples1)

N2 = 32
x_samples2 = np.linspace(0, 1, N2, endpoint=False)
f_samples2 = f(x_samples2)

N3 = 64
x_samples3 = np.linspace(0, 1, N3, endpoint=False)
f_samples3 = f(x_samples3)

N4 = 512
x_samples4 = np.linspace(0, 1, N4, endpoint=False)
f_samples4 = f(x_samples4)

# Plot the sampled function
plt.figure(figsize=(10, 5))
plt.plot(x_samples4, f_samples4, '-', label='f(x)')
plt.plot(x_samples3, f_samples3, '^', label='Sampled f(x) (N=64)')
plt.plot(x_samples2, f_samples2, 's', label='Sampled f(x) (N=32)')
plt.plot(x_samples1, f_samples1, 'o', label='Sampled f(x) (N=16)')
plt.title('Sampled Function f(x)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.grid()
plt.savefig("sampled_function.png")

# Compute the DFT of the sampled function at three values of N
X1 = dft(f_samples1)
k1 = np.arange(N1)

X2 = dft(f_samples2)
k2 = np.arange(N2)

X3 = dft(f_samples3)
k3 = np.arange(N3)

# Sample the analytical Fourier Transform at the same k values for comparison
ksamples = np.linspace(0, N3//2, N4, endpoint=False)
F_samples = F(ksamples)

# Plot the magnitude of the DFT (Only the first N//2 components, the rest are symmetric)
plt.figure(figsize=(10,5))
plt.plot(ksamples, np.abs(F_samples), '-', label='Analytical |F(k)|')
plt.plot(k3[:N3//2], np.abs(X3[:N3//2]), '^', label='N=64')
plt.plot(k2[:N2//2], np.abs(X2[:N2//2]), 's', label='N=32')
plt.plot(k1[:N1//2], np.abs(X1[:N1//2]), 'o', label='N=16')

plt.title('Magnitude of DFT')
plt.xlabel('k')
plt.ylabel('|F[k]|')
plt.grid()
plt.legend()
plt.savefig("dft_magnitude.png")

# Compute the inverse DFT to reconstruct the original function
f_reconstructed1 = idft(X1)
f_reconstructed2 = idft(X2)
f_reconstructed3 = idft(X3)

# Plot the original and reconstructed functions
plt.figure(figsize=(10,5))
plt.plot(x_samples4, f_samples4, '-', label='Original f(x)')
plt.plot(x_samples3, f_reconstructed3.real, '--', label='Reconstructed f(x) (N=64)')
plt.plot(x_samples2, f_reconstructed2.real, '--', label='Reconstructed f(x) (N=32)')
plt.plot(x_samples1, f_reconstructed1.real, '--', label='Reconstructed f(x) (N=16)')

plt.title('Original vs Reconstructed Function')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid()
plt.savefig("original_vs_reconstructed.png")

# Take timing data for the DFT and an off-the-shelf FFT implementation for comparison
import scipy.fft as fft

def time_function(func, data, repeats=50):
    start = time.perf_counter()
    for _ in range(repeats):
        func(data)
    end = time.perf_counter()
    return (end - start) / repeats

dft_times = []
fft_times = []
N_vals = [8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576]

for N in N_vals:
    x_samples = np.linspace(0,1,N,endpoint=False)
    f_samples = f(x_samples)

    if N <= 2048:  # Limit DFT timing to smaller N due to its O(N^2) complexity
        dft_times.append(time_function(dft, f_samples, repeats=3))
    fft_times.append(time_function(fft.fft, f_samples, repeats=2000))

ON2data = [N**2 * dft_times[8]/(N_vals[8] * N_vals[8]) for N in N_vals]
ONlogNdata = [N * np.log2(N) * fft_times[-1]/(N_vals[-1]*np.log2(N_vals[-1])) for N in N_vals]

# Plot the timing data
plt.figure(figsize=(10,5))
plt.plot(N_vals, ON2data, '-', label='O(N^2)') # Plot analytic O(N^2) scaling for DFT
plt.plot(N_vals, ONlogNdata, '--', label=r'O(N $\log_2$ N)') # Plot analytic O(N log N) scaling for FFT
plt.plot([8, 16, 32, 64, 128, 256, 512, 1024, 2048], dft_times, 'o', label='DFT Time')
plt.plot(N_vals, fft_times, 's', label='FFT Time')
plt.title('Timing of DFT vs FFT')
plt.xlabel('N')
plt.ylabel('Time (seconds)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()
plt.savefig("dft_vs_fft_timing.png")