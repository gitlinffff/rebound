# Plot the HST view of dusts with colorcoding of distance to Dimorphos and angle at T0

import os, pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from divide_ejectacone import load_particle_data

def get_r_sky_north():
	# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
	SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
														[ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
														[-0.104839674791979,  0.124915784491013,  0.986612735258626]])

	# calculate sky north vector in 'Sun Body Center' and 
	# convert it to 'Didymos System Barycenter' frame
	obliq_earth = np.deg2rad(23.4392911)
	sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
	r_sky_north = SBC_rotate_DSB @ sky_north
	return r_sky_north

def get_T0_vars(filepath):
	"""
	Calculates dust particles' distance from Dimorphos and their ejection 
	angle relative to the cone axis (negative Y-axis) at T0.
	"""
	data = load_particle_data(filepath)
	p_t0 = data['p_t'][0]

	dimor = p_t0[1]
	dust  = p_t0[4:]
	ids = dust[:, 0].astype(np.int64)

	pos_dust  = dust[:, 1:4]; vel_dust  = dust[:, 4:7]
	pos_dimor = dimor[1:4];   vel_dimor = dimor[4:7]

	# Calculate dust position and velocity relative to Dimorphos
	rel_pos = pos_dust - pos_dimor
	rel_vel = vel_dust - vel_dimor

	# Distance and Velocity Magnitudes
	speed_dsb = np.linalg.norm(vel_dust, axis=1) # relative to DSB
	dists = np.linalg.norm(rel_pos, axis=1)
	speed_dimor = np.linalg.norm(rel_vel, axis=1)# relative to Dimorphos

	# Calculate Theta (Angle relative to -Y axis)
	cone_axis = np.array([0, -1, 0])
	cos_theta = (rel_pos @ cone_axis) / dists
	theta = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

	return {'dists': dists, 'speed_dsb': speed_dsb,
	        'speed_dimor': speed_dimor, 'theta': theta, 'ids': ids}

def process_data(data_dir, RUN_NUMBERS, t0_vars, fraction=1.0):
	"""
	prepare for plotting, combine particle data from all datasets.
	"""
	r_sky_north = get_r_sky_north()

	all_x, all_y, all_dist, all_speed_dimor, all_speed_dsb, all_theta = [], [], [], [], [], []
	for run_idx in RUN_NUMBERS:
		# Load the target dataset (e.g., Day 11)
		filepath = os.path.join(data_dir, f"{run_idx:03d}_snapshots.pkl")
		data = load_particle_data(filepath)
		if data is None:
			return

		particles = data['p_t'][0]
		dusts = particles[4:]
		#rdust = data['radius_dust']
		#day = data['day'][0]

		# randomly select a fraction of all dusts
		Nsample = int(len(dusts) * fraction)
		indices = np.sort(np.random.choice(len(dusts), size=Nsample, replace=False))
		sampled_dusts = dusts[indices, :]
		
		dust_ids = sampled_dusts[:, 0].astype(np.int64)

		# position vector
		r_earth = particles[3, 1:4]  # [x, y, z] of the Earth   (as a proxy of HST)
		r_dust = sampled_dusts[:, 1:4]  # [x, y, z] of all dust particles

		# =============== Project to HST view plane =======================
		# calculate the two basis vectors of the projection plane
		l1 = np.cross(r_sky_north, r_earth)
		l2 = np.cross(r_earth, l1)
		l1 = l1 / np.linalg.norm(l1)
		l2 = l2 / np.linalg.norm(l2)

		# create an array to record coordinates of particles projected onto the plane
		x_proj = r_dust @ l1    # x coordinate
		y_proj = r_dust @ l2    # y coordinate

		# use dust ID to get their initial distance and theta
		mask = np.isin(t0_vars['ids'], dust_ids)
		dist_sel = t0_vars['dists'][mask]
		speed_dimor_sel = t0_vars['speed_dimor'][mask]
		speed_dsb_sel = t0_vars['speed_dsb'][mask]
		theta_sel = t0_vars['theta'][mask]

		all_x.append(x_proj)
		all_y.append(y_proj)
		all_dist.append(dist_sel)
		all_speed_dimor.append(speed_dimor_sel)
		all_speed_dsb.append(speed_dsb_sel)
		all_theta.append(theta_sel)

	return {
			'x': np.hstack(all_x),
			'y': np.hstack(all_y),
			'dist': np.hstack(all_dist),
			'v_d2dimor': np.hstack(all_speed_dimor),
			'v_d2dsb': np.hstack(all_speed_dsb),
			'theta': np.hstack(all_theta)
	}

def plot_t0vars_distribution(data_dict, day_code, output_dir, sort_speed=False):
	day = day_code.split('_')[1]
	alpha = 0.5

	# Accessing dictionary values for brevity
	all_x = data_dict['x']
	all_y = data_dict['y']
	all_dist = data_dict['dist']
	all_v_d2dimor = data_dict['v_d2dimor']
	all_v_d2dsb = data_dict['v_d2dsb']
	all_theta = data_dict['theta']

	if (0):
		print("Plot T0 distance from Dimorphos...", flush=True)
		plt.figure(figsize=(8, 7))
		sca = plt.scatter(all_x, all_y, c=all_dist,
											 s=0.5, cmap='viridis', alpha=alpha,
											vmin=0, vmax=600)
		plt.colorbar(sca, label='Initial Distance from Dimorphos [m]')
		plt.xlim(-2000e3, 2000e3)
		plt.ylim(-2000e3, 2000e3)
		plt.xlabel('HST View Plane X [m]')
		plt.ylabel('HST View Plane Y [m]')
		plt.title(rf'Initial Distance Distribution ($T_0$+{day}day HST FOV)')
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig(os.path.join(output_dir, "init_distance.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
		plt.close()

	if (1):
		print("Plot T0 speed relative to Dimorphos...", flush=True)
		plt.figure(figsize=(8, 7))
		plot_x, plot_y, plot_speed = all_x, all_y, all_v_d2dimor
		if sort_speed:
			# Sort descending (low speed on top of image)
			sort_idx = np.argsort(all_v_d2dimor)[::-1]
			plot_x = all_x[sort_idx]
			plot_y = all_y[sort_idx]
			plot_speed = all_v_d2dimor[sort_idx]
		sca = plt.scatter(plot_x, plot_y, c=plot_speed,
											s=0.5, cmap='tab20', alpha=alpha,
											norm=colors.LogNorm(vmin=0.06, vmax=0.11))
		plt.colorbar(sca, label='Initial speed relative to Dimorphos [m/s]')
		plt.xlim(-2000e3, 2000e3)
		plt.ylim(-2000e3, 2000e3)
		plt.xlabel('HST View Plane X [m]')
		plt.ylabel('HST View Plane Y [m]')
		plt.title(rf'Initial speed Distribution ($T_0$+{day}day HST FOV)')
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig(os.path.join(output_dir, "init_v_d2dimor.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
		plt.close()

	if (1):
		print("Plot T0 speed relative to DSB...", flush=True)
		plt.figure(figsize=(8, 7))
		plot_x, plot_y, plot_speed = all_x, all_y, all_v_d2dsb
		if sort_speed:
			# Sort descending (low speed on top of image)
			sort_idx = np.argsort(all_v_d2dsb)[::-1]
			plot_x = all_x[sort_idx]
			plot_y = all_y[sort_idx]
			plot_speed = all_v_d2dsb[sort_idx]
		sca = plt.scatter(plot_x, plot_y, c=plot_speed,
											s=0.5, cmap='tab20', alpha=alpha,
											norm=colors.LogNorm(vmin=0.15, vmax=0.3))
		plt.colorbar(sca, label='Initial speed relative to DSB [m/s]')
		plt.xlim(-2000e3, 2000e3)
		plt.ylim(-2000e3, 2000e3)
		plt.xlabel('HST View Plane X [m]')
		plt.ylabel('HST View Plane Y [m]')
		plt.title(rf'Initial speed Distribution ($T_0$+{day}day HST FOV)')
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig(os.path.join(output_dir, "init_v_d2dsb.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
		plt.close()

	if (0):
		print("Plot T0 angle within the cone...", flush=True)
		plt.figure(figsize=(8, 7))
		sca = plt.scatter(all_x, all_y, c=all_theta,
											 s=0.5, cmap='viridis', alpha=alpha,
											 vmin=0, vmax=90)
		plt.colorbar(sca, label='Initial Angle in the Cone [degree]')
		plt.xlim(-2000e3, 2000e3)
		plt.ylim(-2000e3, 2000e3)
		plt.xlabel('HST View Plane X [m]')
		plt.ylabel('HST View Plane Y [m]')
		plt.title(rf'Initial Angle Distribution ($T_0$+{day}day HST FOV)')
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig(os.path.join(output_dir, "init_angle.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
		plt.close()
	
	print(f"Images saved to {output_dir}")

if __name__ == "__main__":
	day_code = "day_11.86"
	DATA_DIR = f"/home/linfel/linfel_data/data_high_shortterm_snapshot_data/{day_code}_bldrm"
	RUN_NUMBERS = range(26, 39, 1)
	OUTPUT_DIR = f"/home/linfel/linfel_data/shortterm_anal/{day_code}_bldrm"
	os.makedirs(OUTPUT_DIR, exist_ok=True)

	# get all particles' ID, distance, and theta from T0 snapshot 
	t0_vars = get_T0_vars("/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_0/001_snapshots.pkl")
	# associate particle's coordinates on view plane with their initial distance and theta
	all_data = process_data(DATA_DIR, RUN_NUMBERS, t0_vars, 0.1)
	# create the plot
	plot_t0vars_distribution(all_data, day_code, OUTPUT_DIR, True)
