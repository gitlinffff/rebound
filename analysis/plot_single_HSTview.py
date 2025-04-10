# Plot HST view of particles in several datasets
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from ReadParticle import read_particle_frames
from frame_renderer import init_worker, render_single_HSTview_colorgroups

# data path
data_rootdir = "/u/lli22/ejecta_size_exp"
filenames = [os.path.join(data_rootdir, f"run_{i:03d}", "particles.txt") for i in range(1, 22)]

# Ensure output directory exists for saving frames
output_dir = "/u/lli22/ejecta_size_exp/testplot_Apr4"
os.makedirs(output_dir, exist_ok=True)

datasets = []
# Read particle data
for file in filenames:
    
    Np_seq, time, r_dust, data_p = read_particle_frames(file)
    dataset = {}
    dataset['Np_seq'] = Np_seq
    dataset['time'] = time
    dataset['r_dust'] = r_dust
    dataset['data_p'] = data_p

    datasets.append(dataset)

N_datasets = len(datasets)
print(f"{N_datasets} datasets read in.")


"""HST view of particles in all sizes"""
# parameters
daylist = [131.29]  # set the day to plot
#axis_lims = [-200000e3, 600000e3, -200000e3, 50000e3]    # m
axis_lims = [-50000e3, 50000e3, -40000e3, 30000e3]    # m
size_weight = [1 / N_datasets] * N_datasets
alpha = 0.05  # transparency of the plotted points

# Use a colormap to get distinct colors
colors = plt.cm.tab20.colors + plt.cm.tab20b.colors + plt.cm.tab20c.colors  # rich palette
colors = colors[:N_datasets]  # only take as many as you need
# set up labels
group_labels = {i + 1: [f"{i+1:03d}", '#{:02x}{:02x}{:02x}'.format(
    int(colors[i][0] * 255),
    int(colors[i][1] * 255),
    int(colors[i][2] * 255)
)] for i in range(N_datasets)}

# group_labels = {1: ['1 mm', '#8B4513'],
#                 2: ['3 mm', '#87CEEB'],
#                 3: ['4 mm', '#228B22'],
#                 4: ['5 mm', '#FF4500'],
#                 5: ['7 mm', '#800080'],  # Purple
#                 6: ['1 cm', '#FFD700'],  # Gold (Yellow)
#                 7: ['5 cm', '#00CED1'],
#                 8: ['10 cm', '#DC143C']
#                }
sid = 4  # starting index of dusts (0-3 are Didymos, Dimorphos, Sun, Earth)

# print datetime
start_time = datetime.strptime("2022-09-26 23:17:04.1830", "%Y-%m-%d %H:%M:%S.%f")
new_time = start_time + timedelta(days=daylist[0])
print(new_time)

# check if parameters are correctly set
assert N_datasets==len(size_weight), f"{N_datasets} datasets are given, but the weights of {len(size_weight)} of them are set."

for day in daylist:
    target_seconds = day*86400

    # combine datasets of all sizes
    p_allsize = []
    for i in range(N_datasets):
        # get the frame index at the time
        t_idx = np.searchsorted(datasets[i]['time'], target_seconds, side="left")

        # get the particle data at the time
        p_t = datasets[i]['data_p'][t_idx]

        # add a new column denoting the group of particles
        new_col = np.full((p_t.shape[0], 1), i+1)
        new_col[:sid] = 0  # Didymos, Dimorphos, Sun, Earth should not be grouped with other particles
        p_t = np.hstack((p_t, new_col))

        # randomly select dusts based on weight percentage and add dust data to the list
        num_select = int(size_weight[i] * (len(p_t) - sid))
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
