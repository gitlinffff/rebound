import numpy as np
import matplotlib.pyplot as plt
import miepython as mie

def single_scattering_albedo(r_nm, lambda0, m):


	# Calculate the size parameter (x)
	x = 2 * np.pi * r_nm / lambda0

	# Calculate the efficiencies (Qext, Qsca, etc.)
	qext, qsca, qback, g = mie.efficiencies_mx(m, x)
	
	# Calculate ssa
	ssa = qsca/qext

	return ssa

def plot_qsca(r, qsca, m, lambda0):
# --- Plotting ---
	plt.figure(figsize=(10, 6))

# Plot Qsca vs. radius (r)
	plt.loglog(r, qsca, label=f'm = {m}')

# Set plot labels and title
	plt.xlabel('Particle Radius, $r$ (m)')
	plt.ylabel('Scattering Efficiency, $Q_{\\text{sca}}$')
	plt.title(f'Mie Scattering Efficiency vs. Radius ($\\lambda_0$ = {lambda0:.0f} nm)')
	plt.legend()
	plt.grid(True, which="both", ls="--", linewidth=0.5)

# Optional: Add a vertical line for the wavelength for context
# plt.axvline(lambda0 / (2 * np.pi), color='r', linestyle=':', label='r where x=1')
	plt.savefig("Qsca.png")
#plt.show()

def plot_ssa(r, ssa, m, lambda0):
# --- Plotting ---
	plt.figure(figsize=(10, 6))

# Plot Qsca vs. radius (r)
	plt.plot(np.log10(r), ssa, 'k-', label=f'm = {m}')

# Set plot labels and title
	plt.xlabel('Particle Radius, $r$ (m)')
	plt.ylabel('Single Scattering Albedo $\omega$')
	plt.title(f'($\\lambda_0$ = {lambda0:.0f} nm)')
	plt.legend()
	plt.grid(True, which="both", ls="--", linewidth=0.5)

	plt.savefig(f"ssa_{lambda0:.0f}nm.png")
#plt.show()


def main():
	# wavelength
	wav = np.array([400, 500, 600, 700]) # nm
	
	# refractive index of particle
	m = 1.7 - 0.01j

	# Define the range of particle radii (r)
	radius_min = 1e-6  # m
	radius_max = 1e-1  # m

	# Create a logarithmic array of radii for smooth plotting
	num_points = 200
	r = np.logspace(np.log10(radius_min), np.log10(radius_max), num_points)
	r_nm = r * 1e9  # convert to nm


	for lambda0 in wav:
		ssa = single_scattering_albedo(r_nm, lambda0, m)
		plot_ssa(r, ssa, m, lambda0)



if __name__ == "__main__":
	main()
