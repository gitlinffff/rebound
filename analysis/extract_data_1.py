# Extract snapshot data and save data
import numpy as np
import os, pickle
from ReadParticle import read_specific_frame

def extract(DATA_ROOTDIR, OUTPUT_DIR, TARGET_DAY, RUN_NUMBERS):
	# Calculate target time in seconds
	target_seconds = TARGET_DAY * 86400.

	# Ensure output directory exists for saving frames
	os.makedirs(OUTPUT_DIR, exist_ok=True)
	print(f"Output directory: {OUTPUT_DIR}\n")

	# Loop over the run numbers for consistency
	for run_idx in RUN_NUMBERS:
		input_file = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")

		if not os.path.exists(input_file):
			print(f"Warning: File not found, skipping: {input_file}")
			continue

		Np, time, r_dust, data_p = read_specific_frame(input_file, target_seconds)
	
		# Create a dictionary with the data
		save_data = {
			"day": time / 86400.,
			"Np": Np,
			"radius_dust": r_dust,
			"p_t": data_p
		}

		#output_name = os.path.join(OUTPUT_DIR, f"particle_{run_idx:03d}_day{TARGET_DAY:.2f}.pkl")
		output_name = os.path.join(OUTPUT_DIR, f"{run_idx:03d}_snapshots.pkl")

		# Save the dictionary using pickle
		with open(output_name, 'wb') as f:
				pickle.dump(save_data, f)
		
		print(f"-> Successfully saved to {output_name}\n")
		del Np, time, r_dust, data_p, save_data

def extract_single_day(day):
	# --- Configuration ---
	DATA_ROOTDIR = "/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_run"
	OUTPUT_DIR = f"/home/linfel/linfel_turbo/rebound_exp/data_high_longterm_snapshot_data/day_{day}"
	TARGET_DAY = np.array([float(day)])
	#RUN_NUMBERS = range(33,34) # Define the run numbers you want to process
	RUN_NUMBERS = [25, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40]

	extract(DATA_ROOTDIR, OUTPUT_DIR, TARGET_DAY, RUN_NUMBERS)

def batch_extract_single_day():
	#day_list = ["64.44", "78.65", "83.77", "92.66", "114.75", "131.29", "153.47", "155.31", "177.46", "198.9"]
	day_list = ["78.65", "83.77", "92.66", "114.75", "131.29", "153.47", "155.31", "177.46", "198.9"]
	for day in day_list:
		extract_single_day(day)

if __name__ == "__main__":
	#extract_single_day("198.58")
	batch_extract_single_day()
