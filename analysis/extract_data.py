# Extract data at a time slices and save data
import numpy as np
import os, pickle
from ReadParticle import read_specific_frame

# data path
data_rootdir = "/home/linfel/linfel_scratch/rebound_exp/data_high_shortterm_run_BS"
filenames = [os.path.join(data_rootdir, f"run_{i:03d}", "particles.txt") for i in range(1, 34)]

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/linfel_scratch/rebound_exp/data_high_shortterm_run_BS/day1.18"
os.makedirs(output_dir, exist_ok=True)

# specify time
target_day = 1.18
target_seconds = target_day * 86400.

i = 1
for file in filenames:
	# Read particle data
	Np, time, r_dust, data_p = read_specific_frame(file, target_seconds)

	# Create a dictionary with the data
	save_data = {
		"day": time / 86400.,
		"Np": Np,
		"radius_dust": r_dust,
		"p_t": data_p
	}

	# Create output filename based on dust radius or file index
	output_name = os.path.join(output_dir, f"particle_{i:03d}_day{target_day:.2f}.pkl")

	# Save using pickle
	with open(output_name, 'wb') as f:
		pickle.dump(save_data, f)

	# Optional: clean up memory
	del data_p

	i = i + 1
