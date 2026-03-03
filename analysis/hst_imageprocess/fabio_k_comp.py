import numpy as np
import xarray as xr

def get_compensation_factors(ds, rho=3000, M_kg=3e7, a=2.4):
	"""
	Calculates k_bar_i compensation factors from an xarray Dataset.

	Parameters:
	- ds: xarray.Dataset containing 'radius' and 'particle_count'
	- M_kg: Total estimated mass of the ejecta (default 3e7 kg)
	- rho: Bulk density of particles (default 3500 kg/m3)
	- a: Power law index (default 2.4)
	"""
	# 1. Extract centers and simulated counts
	radii_centers = ds.radius.values
	N_sim_i = ds.particle_count.values

	# 2. Reconstruct Bin Edges from Centers
	# For logarithmic distributions, edges are midpoints in log-space
	log_c = np.log10(radii_centers)
	log_edges = np.zeros(len(log_c) + 1)

	# Calculate midpoints
	log_edges[1:-1] = (log_c[:-1] + log_c[1:]) / 2
	# Extrapolate boundaries
	log_edges[0] = log_c[0] - (log_edges[1] - log_c[0])
	log_edges[-1] = log_c[-1] + (log_c[-1] - log_edges[-2])

	r_binedges = 10**log_edges
	r_m = r_binedges[0]; r_M = r_binedges[-1]

	# 3. Calculate theoretical fraction (ni) per bin
	# ni = integral of r^-a from r_low to r_high
	ni = (r_binedges[1:]**(1-a) - r_binedges[:-1]**(1-a)) / (r_M**(1-a) - r_m**(1-a))

	# 4. Calculate total theoretical particle count (N_bar)
	# Using your derived volumetric integration formula
	sum_term = np.sum(ni * (r_binedges[1:]**2 + r_binedges[:-1]**2) * (r_binedges[1:] + r_binedges[:-1]))
	N_bar = (3 * M_kg) / (rho * np.pi * sum_term)

	# 5. Generate Multiplicative Compensation Factors (k_bar_i)
	if np.any(N_sim_i <= 0):
			print(f"Warning: {np.sum(N_sim_i <= 0)} bins have zero simulated particles.")

	# Concise calculation for k_bar_i
	k_bar_i = np.divide(N_bar * ni, N_sim_i, out=np.zeros_like(ni, dtype=float), where=N_sim_i > 0)

	return k_bar_i, r_binedges

# Usage example:
# k_factors, edges = get_compensation_factors(ds)
# compensated_irradiance = ds.irradiance * xr.DataArray(k_factors, coords={'radius': ds.radius})
