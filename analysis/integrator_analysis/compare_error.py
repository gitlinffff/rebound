import numpy as np
from ReadParticle import read_particle_frames
import matplotlib.pyplot as plt

def get_particle_data(data_path, target_day):
	sid = 4  # index of dusts starting from 4
	target_sec = target_day*86400.
	Np_seq, time, r_dust, data_p = read_particle_frames(data_path)

	t_idx = np.searchsorted(time, target_sec, side="left")
	p_t = data_p[t_idx]

	return p_t[:sid], p_t[sid:] # return mass bodies and particles separately

def find_common_particles(pdata1, pdata2):
	print("finding common particles")
	# Find the particle IDs that are present in BOTH arrays
	common_ids = np.intersect1d(pdata1[:, 0], pdata2[:, 0])

	# Filter both original arrays to keep only the rows with common IDs
	pdata1_filtered = pdata1[np.isin(pdata1[:, 0], common_ids)]
	pdata2_filtered = pdata2[np.isin(pdata2[:, 0], common_ids)]

	return common_ids, pdata1_filtered, pdata2_filtered

def calc_location_l2_error(pdata1, pdata2):
	x_diff = pdata1[:, 1] - pdata2[:, 1]
	y_diff = pdata1[:, 2] - pdata2[:, 2]
	z_diff = pdata1[:, 3] - pdata2[:, 3]
	
	# Calculate the L2 error norm for each particle
	l2_errors = np.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

	return l2_errors, x_diff, y_diff, z_diff

def plot_l2_error_distribution(l2_errors, day, hst_px, output_name):
	# calculate histogram
	bins = np.logspace(np.log10(l2_errors.min()), np.log10(l2_errors.max()), num=100)
	counts, bin_edges = np.histogram(l2_errors, bins=bins)

	# --- Plot histogram ---
	plt.figure(figsize=(8,6))
	plt.bar(bin_edges[:-1], counts, 
					width=bin_edges[1:] - bin_edges[:-1], 
					align='edge', 
					#edgecolor='black',
				 )

	plt.xscale('log')
	plt.xlabel('L2 Error Norm (m)')
	plt.ylabel('Number of Particles')
	plt.title(f'Day {day} Frequency Distribution of Particle Location Errors')
	plt.grid(axis='y', alpha=0.75)

	# Add a vertical line for the mean error and HST pixel size
	plt.axvline(l2_errors.mean(), linestyle='--', color='red', linewidth=2,
							label=f'Mean Error: {l2_errors.mean():.2f} m')
	plt.axvline(hst_px, linestyle='--', color='g', linewidth=2,
							label=f'HST pixel size: {hst_px:.2f} m')
	
	plt.legend()
	plt.savefig(output_name, dpi=150, bbox_inches='tight', pad_inches=0.1)
	#plt.show()

	return

def calc_hst_pixelsize(d):
	"""
	d: distance from HST to target (m)
	return: pixel size in m
	"""
	pixel_arcsec = 0.04
	pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad
	pixel_size = 2 * d * np.tan(pixel_fov/2)

	return pixel_size

if __name__ == '__main__':
	data1_path = "/home/linfel/linfel_scratch/rebound_exp/test_integrator/bs_mm_008/particles.txt"
	data2_path = "/home/linfel/linfel_scratch/rebound_exp/test_integrator/ias_mm_005/particles.txt"
	day = 11.88
	
	_, p1 = get_particle_data(data1_path, day)
	bodies, p2 = get_particle_data(data2_path, day)
	ids, p1, p2 = find_common_particles(p1, p2)
	l2_errors, dx, dy, dz = calc_location_l2_error(p1, p2)
	
	# calculate the scale of one HST pixel based on distance from HST to Didymos
	d = np.linalg.norm(bodies[0,1:4] - bodies[3,1:4])
	px = calc_hst_pixelsize(d)
	
	plot_l2_error_distribution(l2_errors, day, px,
															output_name=f"/home/linfel/linfel_scratch/rebound_exp/test_integrator/postprocess/bs_mm_008cias_mm_005_l2_error_hist_day{day}.png")

	print("\n--- Summary Statistics ---")
	print(f"Mean L2 Error:   {np.mean(l2_errors):.4f}")
	print(f"Median L2 Error: {np.median(l2_errors):.4f}")
	print(f"Max L2 Error:    {np.max(l2_errors):.4f}")
	print(f"Min L2 Error:    {np.min(l2_errors):.4f}")
