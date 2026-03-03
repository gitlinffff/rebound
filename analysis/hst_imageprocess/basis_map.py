import os, pickle, datetime
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from astropy.io import fits
import miepython as mie
from coordinates import r_sky_north, position_dict, day_hstfile_mapping, day_time_mapping
from fabio_k_comp import get_compensation_factors

def process_hst(hst_file, day_code, output_dir, plot_or_not):
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

	if plot_or_not:
		import matplotlib.pyplot as plt
		# Plot HST observation
		plt.figure(figsize=(8, 8))
		pc = plt.pcolormesh(X, Y, log10_hst, cmap='cividis', shading='auto', vmin=-8, vmax=np.nanmax(log10_hst))

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

def pack_simu_results(maps, radii, distances, counts, xedges, yedges, day_code):
	# Calculate centers for the main coordinate axes
	x_centers = (xedges[:-1] + xedges[1:]) / 2
	y_centers = (yedges[:-1] + yedges[1:]) / 2

	ds = xr.Dataset(
		data_vars={
			"irradiance": (["radius", "y", "x"], maps),
			"particle_count": (["radius"], counts),
			"center_distance": (["radius"], distances),
			# Explicitly store the edges
			"xedges": (["x_vertex"], xedges),
			"yedges": (["y_vertex"], yedges),
		},
		coords={
			"radius": radii,
			"x": x_centers,
			"y": y_centers,
		}
	)
	# Optional: Add metadata
	ds.irradiance.attrs['units'] = 'W m-2 um-1 sr-1'
	ds.radius.attrs['units'] = 'm'
	ds.xedges.attrs['units'] = 'm'
	ds.yedges.attrs['units'] = 'm'
	ds.attrs['day_code'] = day_code
	ds.attrs['time_snapshot_str'] = day_time_mapping[day_code]
	ds.attrs['date_created'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

	return ds

def process_simu_intensity(simu_data_dir, RUN_NUMBERS, x_km, y_km, day_code):
	# Parameters
	m = 1.7 - 0.01j       # refractive index of particle
	lambda0 = 555.6e-9    # wavelength in vacuum (m)
	au = 1.495978707e11   # m
	pixel_arcsec = 0.04   # HST pixel size
	pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad

	# Bins of the image
	xedges = x_km * 1e3
	yedges = y_km * 1e3

	n_runs = len(RUN_NUMBERS)
	ny, nx = len(xedges) - 1, len(yedges) - 1
	
	# Pre-allocate NumPy arrays
	# irrad_maps: (n_runs, ny, nx) - optimized for vectorization
	irrad_maps = np.zeros((n_runs, ny, nx))
	d_away = np.zeros(n_runs)
	rlist = np.zeros(n_runs)
	Np_list = np.zeros(n_runs)

	# process every dataset
	for i, run_idx in enumerate(RUN_NUMBERS):
		file = os.path.join(simu_data_dir, f"{run_idx:03d}_snapshots.pkl")
		# Read particle data
		with open(file, 'rb') as f:
			data = pickle.load(f)
			p_t = data['p_t'][0]
			radius_dust = data['radius_dust']
			Np = data['Np'][0]
			day = data['day'][0]
			print(f"using data {file}  day={day}")

		# position vector
		r_earth = p_t[3, 1:4]  # [x, y, z] of the Earth   (as a proxy of HST)
		r_dust = p_t[4:, 1:4]  # [x, y, z] of all dust particles
		r_sun = p_t[2, 1:4]    # [x, y, z] of the Sun

		# =============== particle magnitude =======================
		# vectors and norms
		r_dust2earth = r_earth - r_dust
		r_dust2sun = r_sun - r_dust
		d_PC = np.linalg.norm(r_dust2earth, axis=1)
		d_PS = np.linalg.norm(r_dust2sun, axis=1)

		# Particle single scattering albedo
		size_para = 2 * np.pi * radius_dust / lambda0
		qext, qsca, qback, g = mie.efficiencies_mx(m, size_para)
		omega = qsca/qext

		# phase angle
		cos_alpha = np.einsum('ij,ij->i', r_dust2sun, r_dust2earth) / (d_PC * d_PS) # cos phase angle
		alpha = np.degrees(np.arccos(np.clip(cos_alpha, -1.0, 1.0))) # phase angle (degree)

		# magnitude at 1AU from Sun and observer with 0 phase angle
		m_10 = 5 * np.log10(1329/(2e-3 * radius_dust * np.sqrt(omega)))

		# particle visual magnitude
		vm = m_10 + 5 * np.log10(d_PC*d_PS/au**2) + 0.013*alpha

		# radiometric flux
		E_vega = 3.44e-8  # Vega reference irradiance in visible W m-2 um-1
		E_radio = E_vega * 10 ** ((0.03 - vm)/2.5)  # W m-2 um-1

		# =============== Project to HST view plane =======================
		# calculate the two basis vectors of the projection plane
		l1 = np.cross(r_sky_north, r_earth)
		l2 = np.cross(r_earth, l1)
		l1 = l1 / np.linalg.norm(l1)
		l2 = l2 / np.linalg.norm(l2)

		# create an array to record coordinates of particles projected onto the plane
		x_proj = r_dust @ l1    # x coordinate
		y_proj = r_dust @ l2    # y coordinate

		# ===================== create 2D irradiance map =======================
		# 2D irradiance map
		h_map, _, _ = np.histogram2d(  # irradiance of each pixel (W m-2 um-1)
			x_proj, y_proj,        # x, y coordinates
			bins=[xedges, yedges],
			weights=E_radio,
			density=False  # Set True if you want normalized density
		)
		h_map = h_map.T

		# Consider solid angle of pixel to align the unit as HST data
		# and consider WFC3 F350LP filter
		E_filter = 2.7554e-8  # W m-2 um-1
		irrad_maps[i, :, :] = h_map * (E_filter/E_vega / pixel_fov**2) # W m-2 um-1 sr-1

		# Store scalars
		rlist[i] = radius_dust
		Np_list[i] = Np
		d_away[i] = np.linalg.norm(np.mean(r_dust, axis=0)) # mean distance to DSB
	
	return pack_simu_results(irrad_maps, rlist, d_away, Np_list, xedges, yedges, day_code)		


def synthesize_and_plot(basismap_file, k_factors_path, output_path, use_log=True):
	# Load dataset and apply step slicing (every 10th index)
	with xr.open_dataset(basismap_file) as ds:
		ds_sub = ds.isel(radius=slice(None, None, 1))
		maps = ds_sub.irradiance.values  # (r, y, x)
		radii = ds_sub.radius.values
		x_km, y_km = ds.xedges / 1e3, ds.yedges / 1e3
  
	# Load and map compensation factors
#	with open(k_factors_path, 'rb') as f:
#		k_data = pickle.load(f)
#	bin_idx = np.digitize(radii, k_data['r_binedges']) - 1
#	bin_idx[bin_idx == len(k_data['r_binedges']) - 1] -= 1
#	weights = k_data['k_bar_i'][bin_idx]

	# calculate fabio's compensation factors
	weights, _ = get_compensation_factors(ds_sub, rho=3000, M_kg=3e7, a=2.4)

	# Weighted sum
	synth_map = np.einsum('r,ryx->yx', weights, maps)

	# Toggle Logarithmic Scaling
	if use_log:
		plot_data = np.log10(synth_map + 1e-20) # Add small epsilon to avoid log(0)
		label = r'$\log_{10}$(Brightness) [$W m^{-2} \mu m^{-1} sr^{-1}$]'
		vmax = np.nanmax(plot_data) 
		vmin = vmax - 6  # Show 6 orders of magnitude for best contrast
	else:
		plot_data = synth_map
		label = r'Brightness [$W m^{-2} \mu m^{-1} sr^{-1}$]'
		vmax = np.nanmax(plot_data) / 10.0
		vmin = None

	# Plotting
	fig, ax = plt.subplots(figsize=(8, 8))
	pc = ax.pcolormesh(x_km, y_km, plot_data, cmap='cividis', shading='auto', 
										 vmax=vmax, vmin=vmin)
	cbar = fig.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.5, aspect=30)
	cbar.set_label(label)

	ax.set(xlabel='Projected X [km]', ylabel='Projected Y [km]',
				 title=f'Synthetic HST View ({"Log" if use_log else "Linear"})', aspect='equal')
	ax.grid(True, ls='--', lw=0.1, color='red', alpha=0.7)
	
	plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
	plt.close()

def main(day_code, start, stop):
	# configure file paths
	HST_FILE = os.path.join("/home/linfel/linfel_data/hst_raw_JianyangLi/", day_hstfile_mapping[day_code])

	SIMU_DATA_DIR = ("/home/linfel/linfel_data/"
									 f"data_high_shortterm_snapshot_data/{day_code}_interp")
	RUN_NUMBERS = range(start, stop+1, 1)

	OUTPUT_DIR = f"/home/linfel/linfel_data/shortterm_anal/{day_code}/basis_irrad_maps"
	os.makedirs(OUTPUT_DIR, exist_ok=True)

	# process HST image
	hst_data, log10_hst, x_km, y_km, pixel_km = process_hst(HST_FILE, day_code, OUTPUT_DIR, False)

  # calculate irradiance from simulation results
	ds = process_simu_intensity(SIMU_DATA_DIR, RUN_NUMBERS, x_km, y_km, day_code)
	save_path = os.path.join(OUTPUT_DIR, f"basis_maps_{day_code}_{start}_{stop}.nc")
	ds.to_netcdf(save_path, format="NETCDF4")


if __name__ == "__main__":
  #main("day_11.86", start=10, stop=500)
	synthesize_and_plot("/home/linfel/linfel_data/shortterm_anal/day_11.86/basis_irrad_maps/basis_maps_day_11.86_10_500.nc",
	                    "/home/linfel/linfel_data/shortterm_anal/day_11.86/basis_irrad_maps/k_factors_SNAPSHOT_29.pkl",
											"/home/linfel/linfel_data/shortterm_anal/day_11.86/basis_irrad_maps/day_11.86_compensated_synthetic.png",
											use_log=False)
