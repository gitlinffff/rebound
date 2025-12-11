"""
Use .pkl data
Plot HST view of particles in several datasets
"""

import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from frame_renderer import init_worker, render_single_HSTview_colorgroups

def color_palette(num):
	if num > 60:
		# The 'hsv' colormap provides colors based on the full color wheel (360 degrees)
		cmap = plt.cm.get_cmap('hsv', num)
		cp = [cmap(i)[:3] for i in range(num)]
	else:
		cp = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors) 
		cp = [c[:3] for c in cp[:num]]
	return cp

def plot_HSTview_colorgroups(filenames, output_dir, day, axis_lims):
	"""
	HST view: color groups of particles in different sizes
	"""
	# Set parameters
	N_datasets = len(filenames)
	sid = 4  # starting index of dusts (0-3 are Didymos, Dimorphos, Sun, Earth)

	# Use a colormap to get distinct colors
	colors = color_palette(N_datasets)

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
					p_t = data['p_t']
					radius_dust = data['radius_dust']
					day = data['day']
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


def main_plot():
	day = 5.70  # set the day label
	RUN_NUMBERS = range(200, 271) # select datasets
	
	# data path
	data_dir = (f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/"
							f"day{day:.2f}_interp_run20-27")
	filenames = [os.path.join(data_dir, f"{i:04d}_snapshots.pkl") for i in RUN_NUMBERS]

	# Ensure output directory exists for saving frames
	output_dir = f"/home/linfel/linfel_turbo/rebound_exp/plots/particles_day{day}_interp_data"
	os.makedirs(output_dir, exist_ok=True)

	axis_lims = [-1400e3, 1400e3, -1400e3, 1400e3]    # m
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
	#multi_plots()
	main_plot()
