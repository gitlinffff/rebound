import os, pickle
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
import miepython as mie
import matplotlib.pyplot as plt
from coordinates import position_dict

def process_hst(hst_file, day_code, output_dir):
	# Load FITS image
	with fits.open(hst_file) as hdul:
			hst_data = hdul[0].data  # assuming image is in HDU 0
			header = hdul[0].header

	# Show metadata
	orientat = header.get('ORIENTAT', 'N/A')
	utc_mid = header.get('UTC-MID', 'N/A')
	print(f"HST image metadata:\n orientation: {orientat}\n time: {utc_mid}")

	# clean data, set background to 1e-20, and logarithmic brightness scale
	hst_data[hst_data < 0] = 0
	hst_data[np.isnan(hst_data)] = 0
	log10_hst = np.log10(hst_data + 1e-20)

	# Calculate target distance
	r_Didy_sys_bary = position_dict['Didy_sys_bary'][day_code] # Sun Body Center frame
	r_Hubble = position_dict['Hubble'][day_code]
	target_distance = np.linalg.norm(r_Didy_sys_bary - r_Hubble)  # km

	# Hubble pixel size
	pixel_arcsec = 0.04
	range_km = target_distance
	pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad
	pixel_km = 2 * range_km * np.tan(pixel_fov/2)

	# Axes in km
	ny, nx = hst_data.shape
	x_km = np.arange(nx+1) * pixel_km
	y_km = np.arange(ny+1) * pixel_km
	x_km = x_km - x_km[-1]/2
	y_km = y_km - y_km[-1]/2
	X, Y = np.meshgrid(x_km, y_km)   # km

	# Plot HST observation
	plt.figure(figsize=(8, 8))
	pc = plt.pcolormesh(X, Y, log10_hst, cmap='gray', shading='auto', vmin=-14, vmax=np.nanmax(log10_hst))

	cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.8, aspect=30)  # pad adjusts spacing
	cbar.set_label(r'$\log_{10}$(Pixel Value)')
	
	plt.gca().set_aspect('equal', adjustable='box')
	plt.xlabel('X Distance (km)')
	plt.ylabel('Y Distance (km)')
	plt.title(f'HST Image on {utc_mid}')
	plt.grid(True, linestyle='--', linewidth=0.1, color='red', alpha=0.7)
	plt.tight_layout()

	plt.savefig(os.path.join(output_dir, f'hst_image_{day_code}.png'), dpi=300, bbox_inches='tight')
	#plt.show()
	plt.close()

	return hst_data, log10_hst, x_km, y_km, pixel_km

def process_simu_intensity(simu_data_dir, RUN_NUMBERS, x_km, y_km):
	# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
	SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
														[ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
														[-0.104839674791979,  0.124915784491013,  0.986612735258626]])

	# calculate sky north vector in 'Sun Body Center' and 
	# convert it to 'Didymos System Barycenter' frame
	obliq_earth = np.deg2rad(23.4392911)
	sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
	r_sky_north = SBC_rotate_DSB @ sky_north

	# optical parameters
	m = 1.7 - 0.01j   # refractive index of particle
	lambda0 = 500e-9  # wavelength in vacuum (m)

	xedges = x_km * 1e3
	yedges = y_km * 1e3

	inten_sets = []
	d_away = [] # store distance of cloud of particles to DSB
	rlist = []
	# process every dataset
	for run_idx in RUN_NUMBERS:
		file = os.path.join(simu_data_dir, f"{run_idx:04d}_snapshots.pkl")
		# Read particle data
		with open(file, 'rb') as f:
			data = pickle.load(f)
			p_t = data['p_t']
			radius_dust = data['radius_dust']
			day = data['day']
			print(f"using data {file}  day={day}")

		# =============== Project to HST view plane =======================
		# position vector
		r_earth = p_t[3, 1:4]  # [x, y, z] of the Earth   (as a proxy of HST)
		r_dust = p_t[4:, 1:4]  # [x, y, z] of all dust particles
		#r_sun = p_t[2, 1:4]    # [x, y, z] of the Sun

		# calculate the two basis vectors of the projection plane
		l1 = np.cross(r_sky_north, r_earth)
		l2 = np.cross(r_earth, l1)
		l1 = l1 / np.linalg.norm(l1)
		l2 = l2 / np.linalg.norm(l2)

		# create an array to record coordinates of particles projected onto the plane
		p_projected = np.zeros((len(p_t),3), dtype=float)
		p_projected[:, 0] = p_t[:, 0]                  # copy the column of particle ID
		p_projected[:, 1] = np.dot(p_t[:, 1:4], l1)    # x coordinate
		p_projected[:, 2] = np.dot(p_t[:, 1:4], l2)    # y coordinate

		# rotate particles with an angle on the projected plane
		# This step is to align simulated tail with observation in a line
		#p_projected[:, 1], p_projected[:, 2] = rotate_coords(p_projected[:, 1], p_projected[:, 2], 7.2)
		
		# ===================== Creat 2D density histogram =======================
		hist_data = p_projected[4:]

		# Extract x and y projected coordinates
		x_proj = hist_data[:, 1]
		y_proj = hist_data[:, 2]

		# Mask out particles outside of range
		mask = (
			(x_proj >= xedges[0]) & (x_proj <= xedges[-1]) &
			(y_proj >= yedges[0]) & (y_proj <= yedges[-1])
		)
		x_filtered = x_proj[mask]
		y_filtered = y_proj[mask]

		# Create 2D histogram
		px_den, _, _ = np.histogram2d( # particle density of each pixel
			x_filtered,        # x coordinate
			y_filtered,        # y coordinate
			bins=[xedges, yedges],
			density=False  # Set True if you want normalized density
		)

		# =================== Scattering intensity ========================
		# set weight
		# W = 1e6
		# weight = 10 ** ((np.log10(radius_dust)+1) * np.log10(W) / (-3))    
		weight = 1

		# scattering intensity (assume scattering phase angle constant for all particles)
		x_para = 2*np.pi*radius_dust / lambda0 # size parameter
		qext, qsca, qback, g = mie.efficiencies_mx(m, x_para)
		p_func = 1
		px_inten = qsca * (np.pi * radius_dust**2) * px_den * weight * p_func
		px_inten = px_inten.T

		# ================== Distance of cloud of dust to DSB ====================
		dust_bary = np.mean(r_dust, axis=0)  # average position of the dust particles
		distance_to_DSB = np.linalg.norm(dust_bary)

		# record the results
		inten_sets.append(px_inten)
		rlist.append(radius_dust)
		d_away.append(distance_to_DSB)

	sim_stack = np.stack(inten_sets, axis=-1)            # shape: (ny, nx, n_sizes)
	return sim_stack, np.array(rlist), np.array(d_away)

def residuals(weights, I_models, I_obs):
    # weights can't be negative
    if np.any(weights < 0):
        return np.inf

    # calculate weighted values 
    I_fit = I_models @ weights

    # calculate residuals in log10 space
    epsilon = 1e-20
    log_diff = np.log10(I_obs + epsilon) - np.log10(I_fit + epsilon)
    return log_diff

def fit_weight(simu, obsr):
	# set fitting input X and Y
	N_basis = np.shape(simu)[-1]
	fit_input_X = simu.reshape(-1, N_basis)  # shape: (ny*nx, N_basis)
	fit_input_Y = obsr.flatten()             # shape: (ny*nx)

	# create a mask of HST meaningful signal (pixels of background signal 1e-20 are not used for fitting)
	meaningful_signal_mask = np.log10(fit_input_Y + 1e-20) > -20.0

	# get the meaningful pixels from the HST data and simulated data using the mask
	fit_input_Y_filtered = fit_input_Y[meaningful_signal_mask]
	fit_input_X_filtered = fit_input_X[meaningful_signal_mask, :]

	# --- min/max values ---
	sim_min = np.min(fit_input_X_filtered)
	sim_max = np.max(fit_input_X_filtered)
	hst_min = np.min(fit_input_Y_filtered)
	hst_max = np.max(fit_input_Y_filtered)

	# --- Formatted Print Statement ---
	# The f-string formatting aligns the text and numbers into clean columns.
	# :<20 means left-align in a 20-character space.
	# :>15.2e means right-align in a 15-character space, formatted in scientific notation with 2 decimals.
	print(f"{'':<20} {'min':>15} {'max':>15}")
	print("-" * 55) # Optional: adds a separator line for clarity
	print(f"{'Simulation':<20} {sim_min:>15.2e} {sim_max:>15.2e}")
	print(f"{'HST Observation':<20} {hst_min:>15.2e} {hst_max:>15.2e}")

	# set initial weights
	initial_weights = np.array([1.48607575e-01, 2.03574515e-01, 3.44638220e-01, 4.84672084e-01,
															5.36902768e-01, 5.75193493e-01, 6.08559074e-01, 5.99757568e-01,
															5.61236310e-01, 4.91978292e-01, 3.86263984e-01, 2.45922593e-01,
															7.13588707e-09, 7.13588707e-09])
	initial_weights = np.ones(N_basis)
	
	# least square fitting
	result = least_squares(
			residuals,
			initial_weights,
			args=(fit_input_X_filtered, fit_input_Y_filtered)
	)
	final_weights = result.x

	# --- CALCULATE ERROR BARS (UNCERTAINTIES) ---

	# 1. Extract the Jacobian (J) and the residuals (r)
	J = result.jac
	r = result.fun

	# 2. Define m (data points) and n (parameters)
	m = fit_input_Y_filtered.size
	n = initial_weights.size

	# Degrees of freedom
	dof = m - n

	# 3. Calculate the estimated variance of the residuals (sigma_r^2)
	# The residuals 'r' are in log10 space, so this is the variance of the log_diff.
	r_variance = np.sum(r**2) / dof

	# 4. Calculate the approximated Covariance Matrix (C)
	try:
			J_T_J_inv = np.linalg.inv(J.T @ J)
			C = r_variance * J_T_J_inv
	except np.linalg.LinAlgError:
			# Handle case where the matrix is singular (no unique inverse)
			print("Warning: Cannot calculate covariance matrix. J^T J is singular.")
			weight_errors = np.full_like(final_weights, np.nan)
			I_fit = (fit_input_X @ final_weights).reshape(np.shape(simu)[:2])
			return final_weights, weight_errors, I_fit

	# 5. Extract the standard errors (square root of the diagonal elements)
	weight_errors = np.sqrt(np.diag(C))    

	# calculate weighted fitting values
	I_fit = (fit_input_X @ final_weights).reshape(np.shape(simu)[:2])

	return final_weights, weight_errors, I_fit

def fitting_scatterplot(I_fit, obsr, output_dir):
	x_data = np.log10(I_fit.flatten() + 1e-20)
	y_data = np.log10(obsr.flatten() + 1e-20)

	# Find the combined min and max for both axes
	min_val = min(np.min(x_data), np.min(y_data))
	max_val = max(np.max(x_data), np.max(y_data))

	# Add a little padding to the range
	padding = (max_val - min_val) * 0.05
	axis_min = min_val - padding
	axis_max = max_val + padding
	axis_limits = [axis_min, axis_max]

	bin_edges = np.linspace(axis_min, axis_max, 251)

	# mask 1: Create a mask to find where both I_fit and observation are at their floor value (1e-20)
	zero_mask = (x_data <= -20) | (y_data <= -20.0)
	num_zeros = np.sum(zero_mask)
	x_data_filtered = x_data[~zero_mask]
	y_data_filtered = y_data[~zero_mask]

	# create 2D histogram for performance of model fit
	_den, _, _ = np.histogram2d( # particle density of each pixel
			x_data_filtered,        # x coordinate
			y_data_filtered,        # y coordinate
			bins=[bin_edges, bin_edges],
			density=False  # Set True if you want normalized density
	)
	_den = _den.T

	X, Y = np.meshgrid(bin_edges, bin_edges)

	plt.figure(figsize=(8, 8))
	pc = plt.pcolormesh(X, Y, _den, cmap='cividis', shading='auto')

	# Draw the y=x reference line
	plt.plot(axis_limits, axis_limits, 'r--', label='y=x') # 'r--' is a red dashed line

	# colar bar
	cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.5, aspect=30)
	cbar.set_label(r'Number of Points per Bin')

	plt.gca().set_aspect('equal', adjustable='box')
	plt.xlim(axis_limits)
	plt.ylim(axis_limits)

	plt.xlabel("Predicted Intensity (log10)")
	plt.ylabel("Observed Intensity (log10)")
	plt.title("Model Fit vs. Observation (Density Plot)")
	plt.legend(loc="lower right")
	plt.grid(True, linestyle='--', alpha=0.3)

	plt.savefig(os.path.join(output_dir,'modelfit_scatterplot_.png'), dpi=300, bbox_inches='tight')
	#plt.show()
	plt.close()

def plot_fitted_image(I_fit, pixel_km, output_dir):
	# axes scale
	ny, nx = I_fit.shape
	x_km = np.arange(nx+1) * pixel_km
	y_km = np.arange(ny+1) * pixel_km
	
	# Create meshgrid for bin edges
	X, Y = np.meshgrid(x_km, y_km)   # km

	# Plot pcolormesh
	plt.figure(figsize=(8, 8))

	log10_fit = np.log10(I_fit + 1e-20)
	pc = plt.pcolormesh(X, Y, I_fit, cmap='cividis', shading='auto')

	plt.xlabel('Projected X [km]')
	plt.ylabel('Projected Y [km]')
	plt.title('Intensity Fitting on View Plane')
	plt.grid(True, linestyle='--', linewidth=0.1, color='red', alpha=0.7)
	plt.gca().set_aspect('equal', adjustable='box')

	# Colorbar
	cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.5, aspect=30)
	cbar.set_label(r'$\log_{10}$(Nondimensional Intensity)')

	plt.tight_layout()
	output_name = os.path.join(output_dir, "fitted_image.png")
	plt.savefig(output_name, dpi=300, bbox_inches='tight', pad_inches=0.1)
	#plt.show()
	plt.close()

def plot_w_r(radius, weights, errors, output_dir):
	"""
	plot fitted weights of different dust radius.
	Show error bar for each weight.
	"""
	# calculate the power using a linear regression (A@x=b)
	# assume a power model between weights and radius
	A = np.vstack([np.log10(radius), np.ones(len(radius))]).T
	b = np.log10(weights)
	x = np.linalg.lstsq(A, b, rcond=None)[0] # slope and intercept

	plt.figure(figsize=(8, 6))
	# plot the weights
	plt.loglog(radius, weights, 'bo')
	
	# plot the error bar for weights	
	plt.errorbar(radius, weights, 
							 yerr=errors,
				fmt='bo',    # Blue circles for data points
				#linestyle='-', # Connect the points with a line
				capsize=5,     # Length of the error bar caps
				label=f'Fitted Weights'
			)
	
	# Plot the linear fit line
	fit_line = 10**x[1] * radius**x[0]
	plt.loglog(radius, fit_line, 'r--', label=f'Power Law Fit (slope={x[0]:.2f} intercept={x[1]:.2f})')

	# Add axis labels that reflect the content
	plt.xlabel('Particle Radius (m)')
	plt.ylabel('Fitted Weights')
	plt.title('Weights vs. Radius')
	plt.legend()
	plt.grid(True, which="both", ls="--", linewidth=0.5)
	plt.savefig(os.path.join(output_dir, "w_r.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
	plt.close()

def plot_x_r(radius, distance, output_dir):
	"""
	Plot the relationship between the distance the cloud of particles travel and radius.
	Fit the power law relation.
	"""
	# calculate the power using a linear regression (A@x=b)
	# assume a power model between weights and radius
	A = np.vstack([np.log10(radius), np.ones(len(radius))]).T
	b = np.log10(distance)
	x = np.linalg.lstsq(A, b, rcond=None)[0] # slope and intercept

	plt.figure(figsize=(8, 6))
	# plot the weights
	plt.loglog(radius, distance, 'bo')
	
	# Plot the linear fit line
	fit_line = 10**x[1] * radius**x[0]
	plt.loglog(radius, fit_line, 'r--', label=f'Power Law Fit (slope={x[0]:.2f} intercept={x[1]:.2f})')

	# Add axis labels that reflect the content
	plt.xlabel('Particle Radius (m)')
	plt.ylabel('Distance from Didymos System Barycenter (m)')
	plt.title('Distance traveled vs. Radius')
	plt.legend()
	plt.grid(True, which="both", ls="--", linewidth=0.5)
	plt.savefig(os.path.join(output_dir, "distance_r_1.png"), dpi=150, bbox_inches='tight', pad_inches=0.1)
	plt.close()

def simple_run():
	day_code = "day_5.70"
	hst_file = "/home/linfel/linfel_turbo/hst_raw_JianyangLi/16674/stack_18_long.fits"
	
	simu_data_dir = ("/home/linfel/linfel_turbo/rebound_exp/"
									 "data_high_shortterm_snapshot_data/day5.70_interp_run20-27")
	RUN_NUMBERS = range(200, 271)
	
	output_dir = f"/home/linfel/linfel_turbo/rebound_exp/fit_{day_code}_interp_data_11"
	os.makedirs(output_dir, exist_ok=True)
	
	# process HST image
	hst_data, log10_hst, x_km, y_km, pixel_km = process_hst(hst_file, day_code, output_dir)

	# calculate intensity from simulation results
	sim_stack, radius, distance_away = process_simu_intensity(simu_data_dir, RUN_NUMBERS, x_km, y_km)
	
	# plot distance-away with radius
	#plot_x_r(radius[:-6], distance_away[:-6], output_dir)
	#plot_x_r(radius, distance_away, output_dir)
	#return

	# fit the weights
	weights, errors, I_fit = fit_weight(sim_stack, hst_data)
	
	# output fitting results
	fitting_scatterplot(I_fit, hst_data, output_dir)
	plot_fitted_image(I_fit, pixel_km, output_dir)
	plot_w_r(radius, weights, errors, output_dir)
	np.savetxt(os.path.join(output_dir, "w_r.csv"), np.array([radius, weights, errors]).T, fmt='%.8e', delimiter=',')

def fit_different_regions():
	day_code = "day_64.44"
	hst_file = "/home/linfel/linfel_turbo/hst_raw_JianyangLi/16674/stack_31_long.fits"
	
	output_dir = "/home/linfel/linfel_turbo/rebound_exp/fit_regions_test3"
	os.makedirs(output_dir, exist_ok=True)
	
	# process HST image
	hst_data, log10_hst, x_km, y_km, pixel_km = process_hst(hst_file, day_code, output_dir)

	# calculate intensity from simulation results
	simu_data_dir = "/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data"
	RUN_NUMBERS = range(37, 44)
	sim_stack, radius, distance_away = process_simu_intensity(simu_data_dir, RUN_NUMBERS, x_km, y_km)

	ny, nx = hst_data.shape
	region_dict = {'1': [0, nx, 0, ny],
		             '2': [0, nx, int(0.25*ny), int(0.75*ny)],
	               '3': [0, nx, 0, int(0.75*ny)],
	               '4': [0, nx, int(0.25*ny), ny],
	               '5': [int(0.375*nx), nx, 0, ny],
	               '6': [int(0.375*nx), nx, int(0.25*ny), int(0.75*ny)],
	               '7': [int(0.4*nx), nx, int(0.4*ny), int(0.65*ny)],
								}

	for region_key, region_range in region_dict.items():
		print(f"# working on region {region_key}", flush=True)
		subdir = os.path.join(output_dir, f"region{region_key}")
		os.makedirs(subdir, exist_ok=True)
		
		# crop spatial range
		ix_low = region_range[0]
		ix_up  = region_range[1]
		iy_low = region_range[2]
		iy_up  = region_range[3]
		hst_crop = hst_data[iy_low:iy_up, ix_low:ix_up]
		sim_crop = sim_stack[iy_low:iy_up, ix_low:ix_up, :]
		
		# fit the weights
		weights, errors, I_fit = fit_weight(sim_crop, hst_crop)
		
		# output fitting results
		fitting_scatterplot(I_fit, hst_crop, subdir)
		plot_fitted_image(I_fit, pixel_km, subdir)
		plot_w_r(radius, weights, errors, subdir)
		np.savetxt(os.path.join(subdir, "w_r.csv"), np.array([radius, weights, errors]).T, fmt='%.8e', delimiter=',')

def w_r_from_txt():
	for i in range(1,2):
		print(f"working on No.{i}", flush=True)
		data = np.genfromtxt(f"/home/linfel/linfel_turbo/rebound_exp/fit_regions/region{i}/w_r.csv", delimiter=',')
		plot_w_r(data[:,0], data[:,1], f"/home/linfel/linfel_turbo/rebound_exp/fit_regions/region{i}")

def constrain_mass():
  # simulation data
	simu_data_dir = "/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data"
	RUN_NUMBERS = range(37, 44)

	# density of dust particle
	rho = 3000  # kg/m3
	
	# read fitted weights
	wt_data = np.genfromtxt(f"/home/linfel/linfel_turbo/rebound_exp/fit_regions_test3/region1/w_r.csv", delimiter=',')
	radius_to_weight_map = dict(zip(wt_data[:, 0], wt_data[:, 1]))

	# calculate total mass
	tot_mass = 0.
	for run_idx in RUN_NUMBERS:
		file = os.path.join(simu_data_dir, f"{run_idx:03d}_snapshots.pkl")
		# Read particle data
		with open(file, 'rb') as f:
			data = pickle.load(f)
			Np = data['Np'][0]
			radius_dust = data['radius_dust']
			day = data['day'][0]
		
		wt = radius_to_weight_map[radius_dust]
		print(f"Using data {file}\n day={day:.2f}   Np={Np}   radius={radius_dust:.3e}   weight={wt:.3e}\n")
		tot_mass += Np * wt * rho * 4/3 * np.pi * radius_dust**3
	print(f"Total mass: {tot_mass:.2e} kg")

if __name__ == "__main__":
	simple_run()
	#fit_different_regions()
	#w_r_from_txt()
	#constrain_mass()
