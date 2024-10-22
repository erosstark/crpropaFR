import numpy as np
import matplotlib.pyplot as plt

# Constants for the calculations
Omega_m0 = 0.315  # Matter density parameter
Omega_r0 = 5.373e-5  # Radiation density parameter
n = 1.4  # Model parameter
H0 = 67.4 # Hubble constant in km/s/Mpc

# Calculate R0 (Equation 5)
R0 = - (3 * (3 - n)**2 * H0**2 * Omega_m0) / ((2 * n) * ((n - 3) * Omega_m0 + 2 * (n - 2) * Omega_r0))

# Define the function h(z)
def h(z): # hubbleRate f(R) (Equation 8)
    factor = -2 * n * R0 / (3 * (3 - n)**2 * Omega_m0)
    expression = (n - 3) * Omega_m0 * (1 + z)**(3 / n) + 2 * (n - 2) * Omega_r0 * (1 + z)**((n + 3) / n)
    return (factor * expression) ** 0.5 

print(h(0))
# Generate the x values (z)
z = np.linspace(0, 2.5, 100)

# Calculate the corresponding y values (h(z))
hz = [h(i) for i in z]


def hubbleRate(z):
    omegaM = 0.315
    omegaL = 1 - omegaM
    return H0 * np.sqrt(omegaM * (1 + z)**3 + omegaL)

# Plot the graph
plt.plot(z, hz, label='$f(R)$')
plt.plot(z, hubbleRate(z), label='$\Lambda CDM$')
plt.legend()
plt.xlabel('z')
plt.ylabel('h(z) (km/s/Mpc)')
plt.title('Plot of h(z) from z = 0 to z = 2.5')
plt.grid(True)
plt.show()
    