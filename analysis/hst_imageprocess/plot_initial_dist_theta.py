# Plot the HST view of dusts with colorcoding of distance to Dimorphos and angle at T0

import os, pickle
import numpy as np
import matplotlib.pyplot as plt
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

def get_theta_dist(filepath):
	"""
	Get particles' distance to Dimorphos and angle with respect to cone axis at T0.
	"""

	data_day0 = load_particle_data(filepath)

	particles = data_day0['p_t'][0]
	dimor = particles[1]
	dust = particles[4:]

	ids = dust[:, 0].astype(np.int64)

	r_dust = particles[4:, 1:4]  # [x, y, z] of all dust particles

	negy = np.array([0,-1,0])

	costheta = r_dust @ negy.T / np.linalg.norm(r_dust, axis=1)
	theta = np.degrees(np.arccos(np.clip(costheta, -1.0, 1.0)))

	# distance from Dimorphos
	dist = np.linalg.norm(dust[:, 1:4] - dimor[1:4], axis=1)

	return dist, theta, ids


def process_data(data_dir, RUN_NUMBERS, dist, theta, ids, fraction=1.0):
	"""
	prepare for plotting, combine particle data from all datasets.
	"""
	r_sky_north = get_r_sky_north()

	all_x, all_y, all_dist, all_theta = [], [], [], []
	for run_idx in RUN_NUMBERS:
		# Load the target dataset (e.g., Day 11)
		filepath = os.path.join(data_dir, f"{run_idx:03d}_snapshots.pkl")
		data = load_particle_data(filepath)
		if data is None:
			return

		particles = data['p_t'][0]
		dusts = particles[4:, :]
		#rdust = data['radius_dust']
		#day = data['day'][0]

		# randomly select a fraction of all dusts
		Nsample = int(len(dusts) * fraction)
		indices = np.random.choice(len(dusts), size=Nsample, replace=False)
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
		mask = np.isin(ids, dust_ids)
		dist_sel = dist[mask]
		theta_sel = theta[mask]

		all_x.append(x_proj)
		all_y.append(y_proj)
		all_dist.append(dist_sel)
		all_theta.append(theta_sel)
	
	all_x = np.hstack(all_x)
	all_y = np.hstack(all_y)
	all_dist = np.hstack(all_dist)
	all_theta = np.hstack(all_theta)

	return all_x, all_y, all_dist, all_theta

def plot_dist_theta_distribution(all_x, all_y, all_dist, all_theta, day_code, output_dir):
	day = day_code.split('_')[1]
	
	# Plot initial distance from Dimorphos
	plt.figure(figsize=(8, 7))
	sca = plt.scatter(all_x, all_y, c=all_dist,
	                   s=0.5, cmap='viridis', alpha=0.05,
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

	# Plot initial angle within the cone
	plt.figure(figsize=(8, 7))
	sca = plt.scatter(all_x, all_y, c=all_theta,
	                   s=0.5, cmap='viridis', alpha=0.05,
										 vmin=40, vmax=90)
	plt.colorbar(sca, label='Initial Angle in the Cone [degree]')
	plt.xlim(-2000e3, 2000e3)
	plt.ylim(-2000e3, 2000e3)
	plt.xlabel('HST View Plane X [m]')
	plt.ylabel('HST View Plane Y [m]')
	plt.title(rf'Initial Angle Distribution ($T_0$+{day}day HST FOV)')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(os.path.join(output_dir, "init_angle.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
	plt.close()

if __name__ == "__main__":
	day_code = "day_11.86"
	DATA_DIR = f"/home/linfel/linfel_data/data_high_shortterm_snapshot_data/{day_code}"
	RUN_NUMBERS = range(26, 39, 1)
	OUTPUT_DIR = f"/home/linfel/linfel_data/shortterm_anal/{day_code}"
	os.makedirs(OUTPUT_DIR, exist_ok=True)

	# get all particles' ID, distance, and theta from T0 snapshot 
	dist, theta, ids = get_theta_dist("/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_0/001_snapshots.pkl")
	# associate particle's coordinates on view plane with their initial distance and theta
	all_x, all_y, initial_dist, initial_theta = process_data(DATA_DIR, RUN_NUMBERS, dist, theta, ids, 0.6)
	# create the plot
	plot_dist_theta_distribution(all_x, all_y, initial_dist, initial_theta, day_code, OUTPUT_DIR)
