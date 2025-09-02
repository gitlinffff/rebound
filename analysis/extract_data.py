# Extract data at a time slices and save data
import numpy as np
import os, pickle
from ReadParticle import read_particle_frames

# data path
data_rootdir = "/u/lli22/ejecta_exp_datahigh"
filenames = [os.path.join(data_rootdir, f"run_{i:03d}", "particles.txt") for i in range(1, 22)]

# specify the time of data to extract
target_day = 1.14
target_seconds = target_day * 86400.

# Ensure output directory exists for saving frames
output_dir = f"/u/lli22/ejecta_exp_datahigh/day{target_day:.2f}_pkl"
os.makedirs(output_dir, exist_ok=True)

# Process the data
i = 1
for file in filenames:
    # Read particle data
    Np_seq, time, radii_dust, data_p = read_particle_frames(file)

    t_idx = np.searchsorted(time, target_seconds, side="left")
    p_t = data_p[t_idx]
    day = time[t_idx] / 86400.

    # Create a dictionary with the data
    save_data = {
        "radii_dust": radii_dust,
        "p_t": p_t,
        "day": day
    }

    # Create output filename based on dust radius or file index
    output_name = os.path.join(output_dir, f"particle_{i:03d}_day{target_day:.2f}.pkl")

    # Save using pickle
    with open(output_name, 'wb') as f:
        pickle.dump(save_data, f)
        print(f"data saved to {output_name}")

    # Optional: clean up memory
    del data_p

    i = i + 1

print(f"{i-1} datasets processed.\nExtract data for day {target_day} completed.")
