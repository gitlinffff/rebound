import os, pickle
import numpy as np
from astropy.io import fits
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
	plt.show()
	plt.close()

	return hst_data, log10_hst, x_km, y_km

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
	# process every dataset
	for run_idx in RUN_NUMBERS:
		file = os.path.join(simu_data_dir, f"{run_idx:03d}_snapshots.pkl")
		# Read particle data
		with open(file, 'rb') as f:
			data = pickle.load(f)
			p_t = data['p_t'][0]
			radius_dust = data['radius_dust']
			day = data['day'][0]
			print(f"using data {file}  day={day}")

		# position vector
		r_sun = p_t[2, 1:4]    # [x, y, z] of the Sun
		r_earth = p_t[3, 1:4]  # [x, y, z] of the Earth   (as a proxy of HST)
		r_dust = p_t[4:, 1:4]  # [x, y, z] of all dust particles

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
		
		"""2D density histogram"""
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

		# set weight
		# W = 1e6
		# weight = 10 ** ((np.log10(radius_dust)+1) * np.log10(W) / (-3))    
		weight = 1

		# scattering intensity (assume scattering phase angle constant for all particles)
		qext, qsca, qback, g = mie.efficiencies(m, 2*radius_dust, lambda0)
		p_func = 1
		px_inten = qsca * (np.pi * radius_dust**2) * px_den * weight * p_func
		px_inten = px_inten.T

		inten_sets.append(px_inten)

	sim_stack = np.stack(inten_sets, axis=-1)            # shape: (ny, nx, n_sizes)
	return sim_stack

def main():
	day_code = "day_64.44"
	output_dir = "/home/linfel/linfel_turbo/rebound_exp/plots"
	hst_file = "/home/linfel/linfel_turbo/hst_raw_JianyangLi/16674/stack_31_long.fits"
	hst_data, log10_hst, x_km, y_km = process_hst(hst_file, day_code, output_dir)


	simu_data_dir = "/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data"
	RUN_NUMBERS = range(30, 44)
	process_simu_intensity(simu_data_dir, RUN_NUMBERS, x_km, y_km)
