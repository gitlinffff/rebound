import numpy as np
import matplotlib.pyplot as plt
import miepython as mie

# --- Optical Parameters ---
# Define the range of particle radii (r)
radius_min = 1e-6  # m
radius_max = 1e-1  # m

# Create a logarithmic array of radii for smooth plotting
num_points = 100
r = np.logspace(np.log10(radius_min), np.log10(radius_max), num_points)
r = r * 1e9  # convert to nm

m = 1.7 - 0.01j  # refractive index of particle
lambda0 = 500 # wavelength in vacuum (nm)

# --- Mie Calculation ---

# Calculate the size parameter (x)
# The size parameter is x = 2 * pi * r / lambda0
x = 2 * np.pi * r / lambda0

# Calculate the efficiencies (Qext, Qsca, etc.)
# Note: mie.efficiencies takes the size parameter 'x' (or 2*radius_dust * pi / lambda0),
# NOT just the radius. The line in your original prompt seems to be misusing the 'miepython'
# function signature. The correct size parameter 'x' is used here.
qext, qsca, qback, g = mie.efficiencies_mx(m, x)
#qext, qsca, qback, g = mie.efficiencies(m, 2*radius_dust, lambda0)


# --- Plotting ---
plt.figure(figsize=(10, 6))

# Plot Qsca vs. radius (r)
plt.loglog(r, qsca, label=f'm = {m}')

# Set plot labels and title
plt.xlabel('Particle Radius, $r$ (nm)')
plt.ylabel('Scattering Efficiency, $Q_{\\text{sca}}$')
plt.title(f'Mie Scattering Efficiency vs. Radius ($\lambda_0$ = {lambda0:.0f} nm)')
plt.legend()
plt.grid(True, which="both", ls="--", linewidth=0.5)

# Optional: Add a vertical line for the wavelength for context
# plt.axvline(lambda0 / (2 * np.pi), color='r', linestyle=':', label='r where x=1')
plt.savefig("Qr.png")
#plt.show()
