# Extract data at a time slices and save data
import numpy as np
import os, pickle
from ReadParticle import read_particle_frames

# data path
data_rootdir = "/u/lli22/ejecta_exp_qpr"
filenames = [os.path.join(data_rootdir, f"run_{i:03d}", "particles.txt") for i in range(18, 22)]

# Ensure output directory exists for saving frames
output_dir = "/u/lli22/ejecta_exp_qpr/csv_data/day131.29"
os.makedirs(output_dir, exist_ok=True)

# specify time
target_day = 131.29
target_seconds = target_day * 86400.

i = 1
for file in filenames:
    # Read particle data
    Np_seq, time, radii_dust, data_p = read_particle_frames(file)

    t_idx = np.searchsorted(time, target_seconds, side="left")
    p_t = data_p[t_idx]
    day = time[t_idx] / 86400.

    # save to csv
    output_name = os.path.join(output_dir, f"particle_{i:03d}_day{target_day:.2f}.csv")
    np.savetxt(output_name, p_t, delimiter=',')

    # Optional: clean up memory
    del data_p

    i = i + 1
