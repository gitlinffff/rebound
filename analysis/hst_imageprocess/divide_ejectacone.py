import os, pickle
import numpy as np
import matplotlib.pyplot as plt

def load_particle_data(filepath):
  """
  Loads particle data (ID, position) and radius from a pickled file.
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


def map_radial_particle_groups(filepath):
	"""
	Assign particles to bins based on distance to Dimorphos at T0 + 0 day
	"""
	
	data_day0 = load_particle_data(filepath)

	particles = data_day0['p_t'][0]
	dimor = particles[1]
	dust = particles[4:]
	print(len(dust))
	# distance from Dimorphos
	ids = dust[:, 0].astype(np.int64)
	dist = np.linalg.norm(dust[:, 1:4] - dimor[1:4], axis=1)

	# define bins
	nbins = 5
	#bins = np.logspace(np.log10(dist.min()), np.log10(dist.max()), num=50)
	bins = np.linspace(0, 4000, nbins+1)

	if (0):
		# histogram
		counts, bin_edges = np.histogram(dist, bins=bins)
		
		plt.bar(bin_edges[:-1], counts, 
						color='steelblue', edgecolor='black', alpha=0.8)
		
		#plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('Distance Bins (Radial from Asteroid)')
		plt.ylabel('Number of Particles')
		plt.title(f'Ejecta Particle Distribution\nRange: {dist.min():.1f} to {dist.max():.1f} m')
		plt.grid(axis='y', linestyle='--', alpha=0.6)
		
		plt.tight_layout()
		plt.show()

	# assign to bins
	bin_assignments = np.digitize(dist, bins, right=False)
	bin_assignments[bin_assignments == nbins+1] = nbins

	# Store IDs for each bin
	ids_by_bin = {}
	for i in range(1, nbins+1):
		ids_by_bin[i] = ids[bin_assignments == i]
		print(f"Bin {i}: {len(ids_by_bin[i])} particles identified.")

	return ids_by_bin

def save_binned_particles(new_filepath, ids_by_bin):
	"""
	Loads a new dataset, filters particles by their original IDs, 
	and saves each bin to a separate .pkl file.
	"""
	# Load the target dataset (e.g., Day 11)
	data = load_particle_data(new_filepath)
	if data is None:
		return
	
	particles = data['p_t'][0]
	radius_dust = data['radius_dust']
	day = data['day'][0]
	
	mainbody = particles[:4]
	dust = particles[4:]
	
	# Get the filename without extension to use in the output
	base_name = os.path.splitext(os.path.basename(new_filepath))[0]
	
	for bin_num, stored_ids in ids_by_bin.items():
		# Use np.isin to find particles in this snapshot that were in this bin at Day 0
		mask = np.isin(dust[:, 0].astype(np.int64), stored_ids)
		binned_dust = dust[mask]
		output_particles = np.vstack([mainbody, binned_dust])

		# Prepare the dictionary for the new pickle file
		bin_data = {
			'day': np.array([day]),
			'Np': np.array([len(output_particles)]),
			'radius_dust': radius_dust,
			'p_t': [output_particles],
		}
		
		# Save to file (e.g., bin_1_day_11.86_031_snapshots.pkl)
		output_filename = f"/home/linfel/linfel_data/data_high_shortterm_snapshot_data/make_bins/bin_{bin_num}_{base_name}.pkl"
		with open(output_filename, 'wb') as f:
			pickle.dump(bin_data, f)
				
		print(f"Saved {len(binned_dust)} particles to {output_filename}")

if __name__ == "__main__":
	ids_by_bin = map_radial_particle_groups("/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_0/001_snapshots.pkl")
	save_binned_particles("/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_11.86/019_snapshots.pkl", ids_by_bin)
