# Extract snapshot data and save data
import numpy as np
import os, pickle
from ReadParticle import read_specific_frame

def extract(DATA_ROOTDIR, OUTPUT_DIR, TARGET_DAY, RUN_NUMBERS, bldID_path=None):
	# Calculate target time in seconds
	target_seconds = TARGET_DAY * 86400.

	# Ensure output directory exists for saving frames
	os.makedirs(OUTPUT_DIR, exist_ok=True)
	print(f"Output directory: {OUTPUT_DIR}\n")

	# Load IDs to remove
	remove_ids = None
	if bldID_path and os.path.exists(bldID_path):
		# Assuming IDs are in the first column (usecols=0)
		remove_ids = np.genfromtxt(bldID_path, delimiter=',', usecols=0)
		print(f"Loaded filter file. IDs to remove: {len(remove_ids)}")
	elif bldID_path:
		print(f"Warning: ID file {bldID_path} not found. Proceeding without filtering.")

	# Loop over the run numbers
	for run_idx in RUN_NUMBERS:
		input_file = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")
		if not os.path.exists(input_file):
			print(f"Warning: File not found, skipping: {input_file}")
			continue

		Np, time, r_dust, data_p = read_specific_frame(input_file, target_seconds)
		
		# --- Particle Removal Logic ---
		if remove_ids is not None:
			p_ids = data_p[0][:, 0].astype(int)
			mask = ~np.isin(p_ids, remove_ids)
			data_p[0] = data_p[0][mask]
			Np = len(data_p[0]) 
			print(f"Run {run_idx}: Filtered out {len(mask) - Np} particles.")

		# Create a dictionary with the data
		save_data = {
			"day": time / 86400.,
			"Np": Np,
			"radius_dust": r_dust,
			"p_t": data_p
		}

		# Save the dictionary using pickle
		output_name = os.path.join(OUTPUT_DIR, f"{run_idx:03d}_snapshots.pkl")
		with open(output_name, 'wb') as f:
				pickle.dump(save_data, f)
		
		print(f"-> Successfully saved to {output_name}\n")
		del Np, time, r_dust, data_p, save_data

def extract_single_day(day):
	# --- Configuration ---
	DATA_ROOTDIR = "/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_run"
	OUTPUT_DIR = f"/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data/day_{day}"
	bldID_path = "/home/linfel/linfel_turbo/rebound_exp/datahigh_boulder_hash.csv"
	
	TARGET_DAY = np.array([float(day)])
	RUN_NUMBERS = [25, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40]

	extract(DATA_ROOTDIR, OUTPUT_DIR, TARGET_DAY, RUN_NUMBERS, bldID_path)

def batch_extract_single_day():
	#day_list = ["64.44", "78.65", "83.77", "92.66", "114.75", "131.29", "153.47", "155.31", "177.46", "198.9"]
	day_list = ["78.65", "83.77", "92.66", "114.75", "131.29", "153.47", "155.31", "177.46", "198.9"]
	for day in day_list:
		extract_single_day(day)

if __name__ == "__main__":
	#extract_single_day("198.58")
	batch_extract_single_day()
