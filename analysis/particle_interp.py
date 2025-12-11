import numpy as np
import pickle
import os

def load_particle_data(filepath):
	"""
	Loads particle data (ID, position) and radius from a pickled file.
	
	Args:
			simu_data_dir (str): Directory containing the simulation data files.
			run_idx (int): Index of the simulation run (used for filename).
			
	Returns:
			tuple: (p_t: array of [ID, x, y, z, ...], radius_dust: float)
	"""
	try:
		with open(filepath, 'rb') as f:
			data = pickle.load(f)
			#p_t = data['p_t'][0]  # p_t is 2D array: [ID, x, y, z, vx, vy, vz]
			#radius_dust = data['radius_dust'] # radius_dust is the single radius value for this set
			#return p_t, radius_dust
			return data
	except FileNotFoundError:
		print(f"Error: File not found at {file_path}")
		return None

def interpolate_particle_position(data1, data2, alpha):
	"""
	Interpolates the position and radius of particles between two datasets.
	
	The interpolation is performed only for particles with IDs present in both datasets.
	r_interp = (1-alpha) * r1 + alpha * r2
	P_interp = (1-alpha) * P1 + alpha * P2
	
	Args:
			p1 (np.ndarray): Particle data from first size bin (ID, x, y, z, ...).
			r1 (float): Radius of particles in p1.
			p2 (np.ndarray): Particle data from second size bin (ID, x, y, z, ...).
			r2 (float): Radius of particles in p2.
			alpha (float): Interpolation coefficient (0.0 <= alpha <= 1.0).
			
	Returns:
			tuple: (
					matched_ids: 1D array of particle IDs, 
					p_interp_pos: 2D array of interpolated [x, y, z] coordinates,
					r_interp: float of the interpolated radius
			)
	"""
	if alpha < 0.0 or alpha > 1.0:
		raise ValueError("Interpolation coefficient 'alpha' must be between 0.0 and 1.0")

	r1 = data1['radius_dust']
	p1 = data1['p_t'][0]
	r2 = data2['radius_dust']
	p2 = data2['p_t'][0]
	day = data1['day'][0]

	# Get the ID column (column 0)
	ids1 = p1[:, 0]
	ids2 = p2[:, 0]

	# Find common IDs present in both arrays
	# Return_indices=True helps efficiently retrieve the corresponding data points later
	matched_ids, idx1, idx2 = np.intersect1d(ids1, ids2, assume_unique=True, return_indices=True)

	if len(matched_ids) == 0:
		print("Warning: No common particles found between the two datasets.")
		return np.array([]), np.array([]), np.nan

	# Extract position data (columns 1, 2, 3 for x, y, z)
	p1_pos = p1[idx1, 1:4]
	p2_pos = p2[idx2, 1:4]

	# Calculate interpolated radius and position (P_interp = alpha * P1 + (1 - alpha) * P2)
	r_interp = (1-alpha) * r1 + alpha * r2
	p_interp_pos = (1-alpha) * p1_pos + alpha * p2_pos

	# save data to dictionary
	save_data = {
		"day": day,
		"Np": len(matched_ids),
		"radius_dust": r_interp,
		"p_t": p_interp_pos
	}
	
	return save_data, r1, r2


def run_batch_interp():
	"""
	interpolation between (run_020, run_021), (run_021, run_022), (run_022, run_023), ...
	"""
	simu_data_dir = ("/home/linfel/linfel_turbo/rebound_exp/"
									 "data_high_shortterm_snapshot_data/day5.70")
	RUN_NUMBERS = list(range(20, 27)) 

	output_dir = (f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/"
								f"day5.70_interp_run{RUN_NUMBERS[0]}-{RUN_NUMBERS[-1]+1}")
	os.makedirs(output_dir, exist_ok=True)

	for run_idx in RUN_NUMBERS:
		filepath1 = os.path.join(simu_data_dir, f"{run_idx:03d}_snapshots.pkl")
		filepath2 = os.path.join(simu_data_dir, f"{run_idx+1:03d}_snapshots.pkl")

		# Load data for two neighboring radius bins
		p_data1 = load_particle_data(filepath1)
		p_data2 = load_particle_data(filepath2)

		if p_data1 is not None and p_data2 is not None:
			for alpha in np.linspace(0.0, 0.9, 10):
				# Perform the interpolation
				data_interp, r_1, r_2 = interpolate_particle_position(p_data1, p_data2, alpha)

				# Print results
				print(f"\n{'='*70}")
				print(f"Radius 1: {r_1:.4e} m, Radius 2: {r_2:.4e} m")
				print(f"Interpolation coefficient alpha: {alpha}")
				print(f"Interpolated Radius: {data_interp['radius_dust']:.4e} m")
				print(f"Number of matched particles: {data_interp['Np']}")

				# Save the dictionary using pickle
				output_name = os.path.join(output_dir, f"{run_idx:03d}{int(alpha*10)}_snapshots.pkl")
				with open(output_name, 'wb') as f:
					pickle.dump(data_interp, f)


def run_one():
	"""
	interpolation between (run_020, run_021), (run_021, run_022), (run_022, run_023), ...
	"""
	simu_data_dir = ("/home/linfel/linfel_turbo/rebound_exp/"
									 "data_high_shortterm_snapshot_data/day5.70")
	RUN_NUMBERS = list(range(26, 27)) 

	output_dir = (f"/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_snapshot_data/"
								f"day5.70_interp_run20-27")
	os.makedirs(output_dir, exist_ok=True)

	for run_idx in RUN_NUMBERS:
		filepath1 = os.path.join(simu_data_dir, f"{run_idx:03d}_snapshots.pkl")
		filepath2 = os.path.join(simu_data_dir, f"{run_idx+1:03d}_snapshots.pkl")

		# Load data for two neighboring radius bins
		p_data1 = load_particle_data(filepath1)
		p_data2 = load_particle_data(filepath2)

		if p_data1 is not None and p_data2 is not None:
			alpha=1.0
			# Perform the interpolation
			data_interp, r_1, r_2 = interpolate_particle_position(p_data1, p_data2, alpha)

			# Print results
			print(f"\n{'='*70}")
			print(f"Radius 1: {r_1:.4e} m, Radius 2: {r_2:.4e} m")
			print(f"Interpolation coefficient alpha: {alpha}")
			print(f"Interpolated Radius: {data_interp['radius_dust']:.4e} m")
			print(f"Number of matched particles: {data_interp['Np']}")

			# Save the dictionary using pickle
			output_name = os.path.join(output_dir, f"{run_idx:03d}{int(alpha*10)}_snapshots.pkl")
			with open(output_name, 'wb') as f:
				pickle.dump(data_interp, f)

if __name__ == "__main__":
	run_batch_interp()
	run_one()
