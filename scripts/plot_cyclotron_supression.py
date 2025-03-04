import numpy as np
import matplotlib.pyplot as plt

def calc_Px(x, T):
    """Compute P_x for given x and T."""
    Px_values = np.sqrt(1/x - 1) * np.exp(-511 / T * (1/x - 1))
    Px_values[x >= 1] = 0  # Avoid invalid sqrt() for x >= 1
    return Px_values

def calc_Pnorm_T(T_values, Nx=1001):
    """Compute the integrated normalization Pnorm for a range of temperatures."""
    x_lo, x_hi = 1e-3, 1
    dlogx = np.log10(x_hi / x_lo) / (Nx - 1)
    
    Pnorm = []
    x_vals = x_lo * 10 ** (np.arange(1, Nx) * dlogx)  # Logarithmic grid
    dx_vals = np.diff(x_vals, prepend=x_lo)  # Compute dx
    
    for T in T_values:
        Px_vals = calc_Px(x_vals, T)  # Compute Px at all x
        Ptot = np.sum(Px_vals * dx_vals)  # Midpoint integration
        Ptot = max(Ptot, 1e-10)  # Prevent very small values
        Pnorm.append(Ptot)
    
    return np.array(Pnorm)

# Define x values for plotting Px
x_vals = np.logspace(-3, 0, 1000)  # Logarithmic x from 10^-3 to 1
T_values = [10, 100, 1000, 5000]  # Example temperatures

# Plot P_x for different temperatures
plt.figure(figsize=(8, 5))
for T in T_values:
    plt.plot(x_vals, calc_Px(x_vals, T), label=f"T = {T} keV")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("x (keV)")
plt.ylabel("P(x)")
plt.title("Cyclotron Emission Function P_x vs x")
plt.legend()
plt.grid(True)
plt.show()

# Compute and plot Pnorm(T)
T_range = np.logspace(0, 4, 50)  # T from 1 keV to 10000 keV
Pnorm_values = calc_Pnorm_T(T_range)

plt.figure(figsize=(8, 5))
plt.plot(T_range, Pnorm_values, marker='o', linestyle='-')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("T (keV)")
plt.ylabel("Pnorm(T)")
plt.title("Integrated Cyclotron Emission Normalization Pnorm vs T")
plt.grid(True)
plt.show()

