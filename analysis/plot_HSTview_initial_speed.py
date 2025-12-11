"""
Use .pkl data
Plot HST view with initial speed colorcoding
"""

import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from datetime import datetime, timedelta
from frame_renderer import init_worker, render_single_HSTview_colorgroups

def m_to_km(x, _):
	# Convert length unit for axis labels
	return f'{x / 1e3:.1f}'

def get_speed(filename):
	with open(filename, 'rb') as f:
		data = pickle.load(f)
		p_t = data['p_t'][0]
		radius_dust = data['radius_dust']
		day = data['day'][0]
		print(f"Using data **{filename.split('/')[-1]}**, day={day:.2f}, radius={radius_dust}")
	
	# calculate particle speed (L2 norm of velocity vector)
	speeds = np.linalg.norm(p_t[:, 4:7], ord=2, axis=1)

	return speeds

def project_to_HSTview(p_data):
	# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
	SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
														[ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
														[-0.104839674791979,  0.124915784491013,  0.986612735258626]])
	# calculate sky north vector in 'Sun Body Center' and 
	# convert it to 'Didymos System Barycenter' frame
	obliq_earth = np.deg2rad(23.4392911)
	sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
	r_sky_north = SBC_rotate_DSB @ sky_north

	# position vector of Earth (as a proxy of HST)
	r_earth = p_data[3, 1:4]

	# calculate the two basis vectors of the projection plane
	l1 = np.cross(r_sky_north, r_earth)
	l2 = np.cross(r_earth, l1)
	l1 = l1 / np.linalg.norm(l1)
	l2 = l2 / np.linalg.norm(l2)

	# create an array to record coordinates of particles projected onto the plane
	p_projected = np.zeros((len(p_data),3), dtype=float)
	p_projected[:, 0] = p_data[:, 0]   # copy the column of particle ID
	p_projected[:, 1] = np.dot(p_data[:, 1:4], l1)    # x coordinate
	p_projected[:, 2] = np.dot(p_data[:, 1:4], l2)    # y coordinate

	return p_projected

def plot_colorcoding(p_proj, speed, axis_lim_list, day):
	# Create a new figure for this frame
	fig, ax = plt.subplots(figsize=(8, 8))

	# Determine axis limit
	if axis_lim_list is None:
			print("No axis limits specified, using default limits.")
	elif isinstance(axis_lim_list, list):
			if len(axis_lim_list)==4:         # set upper, bottom, left, right axis limits
					ax.set_xlim(axis_lim_list[0], axis_lim_list[1])
					ax.set_ylim(axis_lim_list[2], axis_lim_list[3])
					ax.set_aspect('equal', adjustable='box')
	elif isinstance(axis_lim_list, float):
			al = axis_lim_list         # Use single float value for a single frame  
			ax.set_xlim(-al, al)
			ax.set_ylim(-al, al)
			ax.set_aspect('equal', adjustable='box')
	
	# plot Didymos and Dimorphos
	ax.scatter(p_proj[0, 1], p_proj[0, 2], c='red', s=10, zorder=3, label='Didymos')
	ax.scatter(p_proj[1, 1], p_proj[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

	# plot dust particles
	p_scat = ax.scatter(p_proj[4:, 1], p_proj[4:, 2], c=speed[4:], cmap='viridis',
											norm=LogNorm(vmin=1e-3, vmax=40.), s=0.5, alpha=0.3)

	# Add a color bar to show the speed scale
	cbar = fig.colorbar(p_scat, ax=ax)
	cbar.set_label(r'$speed$ (m/s)', fontsize=12)

	# Plot Sun direction relative to Didymos System Barycenter
	sun_x, sun_y = p_proj[2, 1], p_proj[2, 2]
	sun_distance = (sun_x**2 + sun_y**2) ** 0.5
	arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
	arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
	ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')
	
	# Customize axis
	ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
	ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
	ax.set_xlabel('x / km')
	ax.set_ylabel('y / km')
	ax.grid()
	ax.legend(loc='upper right')
 
	# Set title
	ax.set_title(f't = T0 + {day:.2f} days   HST Perspective')

	return fig

def main():
	"""
	Plot particles with initial speed colorcoding.
	Input:
		speed_data: 1d array of speed of particles
	"""
	speedfile = f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/day0/001_snapshots.pkl"
	#pfile = f"/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data/035_snapshots.pkl"
	pfile = speedfile

	output_dir = f"/home/linfel/linfel_turbo/rebound_exp/plots/plots_day64.44"
	os.makedirs(output_dir, exist_ok=True)

	axis_lims = [-10000e3, 10000e3, -10000e3, 10000e3]    # m
	axis_lims = 3e3    # m
	
	with open(pfile, 'rb') as f:
		data = pickle.load(f)
		p_t = data['p_t'][0]
		radius_dust = data['radius_dust']
		day = data['day'][0]
		print(f"Using data **{pfile.split('/')[-1]}**, day={day:.2f}, radius={radius_dust}")

	# get particles' initial speeds
	speed_data = get_speed(speedfile)
	
	# match the particles with their initial speed using their IDs.
	p_id = p_t[:, 0]
	idx = p_id.astype(int) - 1
	speed_sel = speed_data[idx]
	
	# project to HST view 
	projected_particles = project_to_HSTview(p_t)
	
	# make the plot
	fig = plot_colorcoding(projected_particles, speed_sel, axis_lims, day)

	# Save the image
	frame_filename = os.path.join(output_dir, f"hstview_init_speed_day{day:.2f}.png")
	fig.savefig(frame_filename, dpi=300, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig)


#======================================================================================


def plot_HSTview_colorgroups(filenames, output_dir, day, axis_lims):
	"""
	HST view: color groups of particles in different sizes
	"""
	# Set parameters
	N_datasets = len(filenames)
	sid = 4  # starting index of dusts (0-3 are Didymos, Dimorphos, Sun, Earth)

	# Use a colormap to get distinct colors
	colors = plt.cm.tab20.colors + plt.cm.tab20b.colors + plt.cm.tab20c.colors  # rich palette
	colors = colors[:N_datasets]  # only take as many as you need

	# print datetime
	start_time = datetime.strptime("2022-09-26 23:17:04.1830", "%Y-%m-%d %H:%M:%S.%f")
	new_time = start_time + timedelta(days=day)
	print(f"Plotting T0+{day} day, {new_time}")
	target_seconds = day * 86400.

	# combine datasets of all sizes
	p_allsize = []
	group_labels = {}

	for i in range(N_datasets):
			# Read particle data
			with open(filenames[i], 'rb') as f:
					data = pickle.load(f)
					p_t = data['p_t'][0]
					radius_dust = data['radius_dust']
					day = data['day'][0]
					print(f"Using data **{filenames[i].split('/')[-1]}**, day={day}, radius={radius_dust}")

			# set group label
			group_labels[i+1] = [f"{radius_dust:.2e}", '#{:02x}{:02x}{:02x}'.format(
					int(colors[i][0] * 255),
					int(colors[i][1] * 255),
					int(colors[i][2] * 255))]

			# add a new column denoting the group of particles
			new_col = np.full((p_t.shape[0], 1), i+1)
			new_col[:sid] = 0  # Didymos, Dimorphos, Sun, Earth should not be grouped with other particles
			p_t = np.hstack((p_t, new_col))

			# add Didymos, Dimorphos, Sun, Earth to the beginning of the list
			if i==0:  #either average the position of sun and earth across datasets, or interpolate position of dusts to a certain moment
					p_allsize.append(p_t[:sid])

			# randomly select dusts based on weight percentage and add dust data to the list
			size_weight = 1 / N_datasets  # custom weight
			num_select = int(size_weight * (len(p_t) - sid))
			selected_indices = np.random.choice(range(sid,len(p_t)), num_select, replace=False)
			p_allsize.append(p_t[selected_indices])

	p_allsize = np.vstack(p_allsize)

	# pass global variables
	init_worker(data_p=[p_allsize], time=np.array([target_seconds]), axis_lim_list=axis_lims, output_dir=output_dir, group_labels=group_labels)
	# render the image
	alpha = 0.05 # transparency of the plotted points
	render_single_HSTview_colorgroups(0, alpha=alpha)


def main_plot(RUN_NUMBERS):
	day = 14.91  # set the day label
	RUN_NUMBERS = range(30,44) # select datasets
	
	# data path
	data_dir = f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/day{day}"
	filenames = [os.path.join(data_dir, f"{i:03d}_snapshots.pkl") for i in RUN_NUMBERS]

	# Ensure output directory exists for saving frames
	output_dir = f"/home/linfel/linfel_turbo/rebound_exp/plots/particles_day{day}"
	os.makedirs(output_dir, exist_ok=True)

	axis_lims = [-2000e3, 2000e3, -2000e3, 2000e3]    # m
	#axis_lims = [-28799e3, 28799e3, -28799e3, 28799e3]    # m
	#axis_lims = [-600e3, 600e3, -600e3, 600e3]

	plot_HSTview_colorgroups(filenames, output_dir, day, axis_lims)

def multi_plots():
	day = 5.70  # set the day label
	RUN_RANGES = [   # select data run ranges
        #(25, 42),  # Runs 25 through 41
        #(22, 40),
				(20, 36)
				]
	
	for start_num, stop_num in RUN_RANGES:
		# data path
		data_dir = f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/day{day:.2f}"
		filenames = [os.path.join(data_dir, f"{i:03d}_snapshots.pkl") for i in range(start_num, stop_num)]

		# Ensure output directory exists for saving frames
		output_dir = (f"/home/linfel/linfel_turbo/rebound_exp/plots/"
									f"particles_day{day}_runs{start_num}-{stop_num-1}")
		os.makedirs(output_dir, exist_ok=True)

		axis_lims = [-1400e3, 1400e3, -1400e3, 1400e3]    # m
		#axis_lims = [-2000e3, 2000e3, -2000e3, 2000e3]    # m
		#axis_lims = [-28799e3, 28799e3, -28799e3, 28799e3]    # m
		#axis_lims = [-600e3, 600e3, -600e3, 600e3]

		plot_HSTview_colorgroups(filenames, output_dir, day, axis_lims)

if __name__ == "__main__":
	main()
