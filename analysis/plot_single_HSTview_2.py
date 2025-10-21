"""
Use .pkl data
Plot HST view of particles in several datasets
"""

import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from ReadParticle import read_particle_frames
from frame_renderer import init_worker, render_single_HSTview_colorgroups

# data path
data_rootdir = "/home/linfel/linfel_turbo/rebound_exp/snapshot_data/day11.88"
filenames = [os.path.join(data_rootdir, f"particle_{i:03d}_day11.88.pkl") for i in [1,11,21,31,40,50]]

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/linfel_turbo/rebound_exp/plots/plots_day11.88"
os.makedirs(output_dir, exist_ok=True)

"""HST view of particles in all sizes"""
# Set parameters
day = 11.88                # set the day to plot
#axis_lims = [-10000e3, 35000e3, -5000e3, 2000e3]    # m
#axis_lims = [-28799e3, 28799e3, -28799e3, 28799e3]    # m
axis_lims = [-600e3, 600e3, -600e3, 600e3]
alpha = 0.05                # transparency of the plotted points
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
        p_t = data['p_t']
        radius_dust = data['radius_dust']
        day = data['day']
        
    # set group label
    group_labels[i+1] = [f"{radius_dust:.2e}", '#{:02x}{:02x}{:02x}'.format(
        int(colors[i][0] * 255),
        int(colors[i][1] * 255),
        int(colors[i][2] * 255))]

    # add a new column denoting the group of particles
    new_col = np.full((p_t.shape[0], 1), i+1)
    new_col[:sid] = 0  # Didymos, Dimorphos, Sun, Earth should not be grouped with other particles
    p_t = np.hstack((p_t, new_col))

    # randomly select dusts based on weight percentage and add dust data to the list
    size_weight = 1 / N_datasets  # custom weight
    num_select = int(size_weight * (len(p_t) - sid))
    selected_indices = np.random.choice(range(sid,len(p_t)), num_select, replace=False)
    p_allsize.append(p_t[selected_indices])

    # add Didymos, Dimorphos, Sun, Earth to the beginning of the list
    if i==0:  #either average the position of sun and earth across datasets, or interpolate position of dusts to a certain moment
        p_allsize.insert(0, p_t[:sid])

p_allsize = np.vstack(p_allsize)

# pass global variables
init_worker(data_p=[p_allsize], time=np.array([target_seconds]), axis_lim_list=axis_lims, output_dir=output_dir, group_labels=group_labels)
# make plot for single frame
render_single_HSTview_colorgroups(0, alpha=alpha)
