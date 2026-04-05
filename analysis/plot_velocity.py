import os, pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def load_particle_data(filepath):
  """
  Loads particle data (ID, position) and radius from a pickled file.
  """
  try:
    with open(filepath, 'rb') as f:
      data = pickle.load(f); return data
  except FileNotFoundError:
    print(f"Error: File not found at {file_path}"); return None

def decompose_velocity(v_vecs, n_dir):
	"""decompose to parallel and perpendicular components"""
	n_hat = n_dir / np.linalg.norm(n_dir)
	v_para = np.dot(v_vecs, n_hat)
	v_perp_vec = v_vecs - np.outer(v_para, n_hat)
	v_perp = np.linalg.norm(v_perp_vec, axis=1)	
	return np.column_stack((v_para, v_perp))

def get_T0_vars(filepath):
	data = load_particle_data(filepath)
	p_t = data['p_t'][0]
	
	# Extract states
	dimor_v = p_t[1][4:7]
	sun_pos = p_t[2][1:4]
	dust_data = p_t[4:]
	
	# Directions
	n_srp = -sun_pos / np.linalg.norm(sun_pos)
	n_imp = np.array([0, 1, 0])
	
	# Raw Velocities
	v_dsb = dust_data[:, 4:7]
	v_dimor = v_dsb - dimor_v
	
	# Decompositions
	dist_dimor_imp = decompose_velocity(v_dimor, n_imp)
	dist_dimor_srp = decompose_velocity(v_dimor, n_srp)
	dist_dsb_imp   = decompose_velocity(v_dsb, n_imp)
	dist_dsb_srp   = decompose_velocity(v_dsb, n_srp)
	
	# Create DataFrame: using dict to handle 2D arrays (N, 2) properly
	df = pd.DataFrame({
			'ids': dust_data[:, 0].astype(np.int64),
			'v_dsb': list(v_dsb),
			'v_dimor': list(v_dimor),
			'speed_dsb': np.linalg.norm(v_dsb, axis=1),
			'speed_dimor': np.linalg.norm(v_dimor, axis=1),
			'v_dist_dimor_imp': list(dist_dimor_imp),
			'v_dist_dimor_srp': list(dist_dimor_srp),
			'v_dist_dsb_imp': list(dist_dsb_imp),
			'v_dist_dsb_srp': list(dist_dsb_srp)
	})
	return df.set_index('ids')

def process_data(data_dir, run_idx, t0_df):
	filepath = os.path.join(data_dir, f"{run_idx:03d}_snapshots.pkl")
	data = load_particle_data(filepath)
	if not data: return None

	current_ids = data['p_t'][0][4:, 0].astype(np.int64)
	filtered_df = t0_df.loc[t0_df.index.intersection(current_ids)]
	
	return filtered_df

def process_and_export_csv(data_dir, run_idx, t0_vars, output_path):
	filepath = os.path.join(data_dir, f"{run_idx:03d}_snapshots.pkl")
	data = load_particle_data(filepath)

	# 1. Filter particles present in current snapshot
	current_ids = data['p_t'][0][4:, 0].astype(np.int64)
	mask = np.isin(t0_vars['ids'], current_ids)

	# 2. Extract state vectors
	v_dsb = t0_vars['v_dsb'][mask]
	v_dimor = t0_vars['v_dimor'][mask]
	s_dsb = t0_vars['speed_dsb'][mask]
	s_dimor = t0_vars['speed_dimor'][mask]

	# 3. Create a DataFrame for ParaView
	df = pd.DataFrame({
			'vx_dsb': v_dsb[:, 0], 'vy_dsb': v_dsb[:, 1], 'vz_dsb': v_dsb[:, 2],
			'vx_dimor': v_dimor[:, 0], 'vy_dimor': v_dimor[:, 1], 'vz_dimor': v_dimor[:, 2],
			'speed_dsb': s_dsb,
			'speed_dimor': s_dimor,
			# Threshold logic (1 for True, 0 for False)
			'dsb_gt_025': (s_dsb > 0.25).astype(int),
			'dimor_gt_015': (s_dimor > 0.15).astype(int)
	})

	df.to_csv(output_path, index=False)
	print(f"Exported {len(df)} particles to {output_path}")

def plot_velocity_space(vel_vecs, speeds, title, out_path):
	fig = plt.figure(figsize=(10, 8))
	ax = fig.add_subplot(111, projection='3d')
	
	sc = ax.scatter(vel_vecs[:, 0], vel_vecs[:, 1], vel_vecs[:, 2], 
									c=speeds, cmap='viridis', s=1, alpha=0.6,
									norm=colors.LogNorm(vmin=0.01, vmax=20))
	
	ax.set_xlabel('$v_x$ (m/s)'); ax.set_ylabel('$v_y$ (m/s)'); ax.set_zlabel('$v_z$ (m/s)')
	ax.set(xlim=(-20, 20), ylim=(-20, 20), zlim=(-20, 20))
	ax.set_title(title)
	plt.colorbar(sc, label='Speed (m/s)')
	
	plt.savefig(out_path, dpi=300, bbox_inches='tight')
	plt.show()
	plt.close()

def plot_velocity_space_2(vel_vecs, speeds, s_bound, title, out_path):
	fig = plt.figure(figsize=(10, 8))
	ax = fig.add_subplot(111, projection='3d')
	
	low_speed = speeds < s_bound
	ax.scatter(*vel_vecs[~low_speed].T, color='steelblue', s=1, alpha=0.2,
	           label=rf'$v \geq {s_bound}$ m/s')
	ax.scatter(*vel_vecs[low_speed].T, color='darkorange', s=3, alpha=0.9,
	          label=rf'$v < {s_bound}$ m/s')

	ax.set(xlim=(-20, 20), ylim=(-20, 20), zlim=(-20, 20), 
				 xlabel='$v_x$ (m/s)', ylabel='$v_y$ (m/s)', zlabel='$v_z$ (m/s)', title=title)
	ax.legend(loc='upper right', markerscale=8)
	plt.savefig(out_path, dpi=300, bbox_inches='tight')
	plt.close()

def plot_velocity_projections(vel_vecs, speeds, s_bound, title, out_path):
	fig, axes = plt.subplots(1, 3, figsize=(16, 5))

	lim = 3
	low_speed = speeds < s_bound
	# Pairs for vx-vy, vy-vz, vx-vz
	planes = [(0, 1, '$v_x$', '$v_y$'), (1, 2, '$v_y$', '$v_z$'), (0, 2, '$v_x$', '$v_z$')]
	
	for ax, (idx1, idx2, lbl1, lbl2) in zip(axes, planes):
		# Plot high speed
		ax.scatter(vel_vecs[~low_speed, idx1], vel_vecs[~low_speed, idx2], 
							 color='steelblue', s=1, alpha=0.1, label=rf'$v \geq {s_bound}$ m/s')
		# Plot low speed
		ax.scatter(vel_vecs[low_speed, idx1], vel_vecs[low_speed, idx2], 
							 color='darkorange', s=3, alpha=0.6, label=rf'$v < {s_bound}$ m/s')

		ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel=f'{lbl1} (m/s)', ylabel=f'{lbl2} (m/s)')
		ax.set_aspect('equal', adjustable='box')
		ax.grid(True, linestyle='--', alpha=0.5)
		ax.legend(loc='upper right', markerscale=5)

	fig.suptitle(title)
	plt.tight_layout()
	plt.savefig(out_path, dpi=300)
	plt.close()

def plot_para_perp_dist(v_para, v_perp, speeds, s_bound, title, out_path, xlabel):
	plt.figure(figsize=(8, 5))
	lim = 3
	low_speed = speeds < s_bound

	plt.scatter(v_para[~low_speed], v_perp[~low_speed], color='steelblue',
	            s=1, alpha=0.2, label=rf'$v \geq {s_bound}$ m/s')
	plt.scatter(v_para[low_speed], v_perp[low_speed], color='darkorange',
	            s=3, alpha=0.8, label=rf'$v < {s_bound}$ m/s')

	plt.axvline(0, color='black', lw=1, ls='--')
	plt.xlabel(f'{xlabel} (m/s)')
	plt.ylabel('$v_{\perp}$ (m/s)')
	plt.legend(loc='upper right', markerscale=5)
	plt.title(title)
	plt.gca().set_aspect('equal', adjustable='box')
	plt.xlim(-lim, lim)
	plt.ylim(0, lim)
	plt.grid(True, alpha=0.3)
	plt.savefig(out_path, dpi=300)
	plt.close()

def plot_para_perp_density(v_para, v_perp, speeds, s_bound, title, out_path, xlabel):
	plt.figure(figsize=(8, 5))
	lim = 0.5
	resolution = 0.001 # m/s per bin
	nx = int((2 * lim) / resolution); ny = int(lim / resolution)
	
	# Create the 2D Histogram (Density Map)
	plt.hist2d(v_para, v_perp, bins=[nx,ny], range=[[-lim, lim], [0, lim]], 
						 cmap='viridis', norm=colors.LogNorm())
	
	# Overlay the speed boundary as a visual guide (a quarter circle)
	theta = np.linspace(0, np.pi, 100)
	plt.plot(s_bound * np.cos(theta), s_bound * np.sin(theta), 
					 color='red', linestyle='--', lw=1, label=f'Threshold: {s_bound} m/s')

	plt.colorbar(label='Particle Density (log scale)')
	plt.axvline(0, color='white', lw=0.8, ls='--', alpha=0.5)
	plt.xlabel(f'{xlabel} (m/s)')
	plt.ylabel('$v_{\perp}$ (m/s)')
	plt.title(title)
	plt.legend(loc='upper right')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.grid(True, alpha=0.2, color='white')
	plt.savefig(out_path, dpi=300, bbox_inches='tight')
	plt.close()

# --- Execution ---
OUTPUT_DIR = "/home/linfel/linfel_data/shortterm_anal/day_11.86_bldrm"
DATA_DIR = "/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_11.86_bldrm"
T0_FILE = "/home/linfel/linfel_data/data_high_shortterm_snapshot_data/day_0/001_snapshots.pkl"
os.makedirs(OUTPUT_DIR, exist_ok=True)

t0_df = get_T0_vars(T0_FILE)
df_filtered = process_data(DATA_DIR, 33, t0_df)

# --- Configuration Set ---
plot_configs = {  # dust velocity is relative to Dimor or DSB
	'Dimorphos': {
		'v_vecs': np.stack(df_filtered['v_dimor'].values),
		'speed': df_filtered['speed_dimor'].values,
		'v_dist_imp': np.stack(df_filtered['v_dist_dimor_imp'].values),
		'v_dist_srp': np.stack(df_filtered['v_dist_dimor_srp'].values),
		's_bound': 0.15,
		'prefix': 'dimor'
	},
	'DSB': {
		'v_vecs': np.stack(df_filtered['v_dsb'].values),
		'speed': df_filtered['speed_dsb'].values,
		'v_dist_imp': np.stack(df_filtered['v_dist_dsb_imp'].values),
		'v_dist_srp': np.stack(df_filtered['v_dist_dsb_srp'].values),
		's_bound': 0.2504,
		'prefix': 'dsb'
	}
}

# --- Execution Loop ---
for frame, cfg in plot_configs.items():
	# 1. 1x3 Projections
	plot_velocity_projections(
		cfg['v_vecs'], cfg['speed'], cfg['s_bound'],
		f"{frame}-Relative Velocity Projections",
		os.path.join(OUTPUT_DIR, f"{cfg['prefix']}_v_projections.png")
	)

	# 2. Parallel/Perpendicular relative to Impact
	plot_para_perp_density(
		cfg['v_dist_imp'][:, 0], cfg['v_dist_imp'][:, 1], 
		cfg['speed'], cfg['s_bound'],
		f"Ejecta Velocity Distribution (Impact Direction) {frame}-rel.",
		os.path.join(OUTPUT_DIR, f"{cfg['prefix']}_v_dist_impact.png"),
		"$v_{\parallel, \t{Impact}}$"
	)

	# 3. Parallel/Perpendicular relative to SRP
	plot_para_perp_density(
		cfg['v_dist_srp'][:, 0], cfg['v_dist_srp'][:, 1], 
		cfg['speed'], cfg['s_bound'],
		f"Ejecta Velocity Distribution (SRP Direction) {frame}-rel.",
		os.path.join(OUTPUT_DIR, f"{cfg['prefix']}_v_dist_srp.png"),
		"$v_{\parallel, \t{SRP}}$"
	)

print(f"Success: All plots saved to {OUTPUT_DIR}")
