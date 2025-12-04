"""
Quick plot top view, HST view and side view
"""

import os, glob
import pickle
import numpy as np
from ReadParticle import read_specific_frame, read_meta
from frame_renderer import init_worker, render_frame_topview, render_frame_HSTview, render_frame_sideview

#------------------------------------------------------------------------------
def quickview_3(data_path, output_dir = '.', target_days=None, axis_lim=3e3):
	# Read particle data
	target_seconds = np.array(target_days) * 86400.
	Np, time, r_dust, data_p = read_specific_frame(data_path, target_seconds)
	
	# Ensure output directory exists for saving frames
	os.makedirs(output_dir, exist_ok=True)
	
	for i in range(len(time)):
		init_worker(data_p, time, axis_lim[i], output_dir)
		
		render_frame_topview(i)
		render_frame_sideview(i)
		render_frame_HSTview(i)

def quick_plot_particles():
	DATA_ROOTDIR = "/home/linfel/linfel_scratch/rebound_exp/data_high_longterm_run_034-050"
	RUN_NUMBERS = range(34, 51) # Define the run numbers you want to process
	PLOT_DAYS = [0, 1.14, 5.7, 8.8, 11.86, 14.9, 131.29, 153.47]
	AXIS_LIM = [3e3, 200e3, 200e3, 200e3, 300e3, 300e3, 20000e3, 25000e3]        # correspond to each day

	for run_idx in RUN_NUMBERS:
		print(f"processing run_{run_idx:03d}")
		data_path = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")
		output_dir = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "quickview")
		quickview_3(data_path, output_dir, PLOT_DAYS, AXIS_LIM)

def check_metadata():
	"""
	print out time in days of data particles.txt
	"""
	DATA_ROOTDIR = "/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_run_BS"
	RUN_NUMBERS = range(20,25) # Define the run numbers you want to process
	
	for run_idx in RUN_NUMBERS:
		print(f"{'='*70}\nProcessing run_{run_idx:03d}")
		data_path = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")
		Np, time = read_meta(data_path)
		print(f"# of snapshots: {len(time)}")
		print(f"First 5 time: {time[:5]/86400.} day")
		print(f"Last  5 time: {time[-5:]/86400.} day")

#------------------------------------------------------------------------------
def check_pkl_structure(folder_path):
	"""
	Finds all .pkl files in the specified folder and processes them
	to check the structure of 'Np', 'p_t', 'radius_dust', and 'day' keys.
	"""
	pkl_files = glob.glob(os.path.join(folder_path, '*.pkl'))
	pkl_files = sorted(pkl_files)

	if not pkl_files:
		print(f"No .pkl files found in '{folder_path}'.")
		return
	
	for datapath in pkl_files:
		print(f"\n{'='*70}")
		print(f"Processing file: **{os.path.basename(datapath)}**")
		print(f"{'='*70}")

		try:
			with open(datapath, 'rb') as f:
				data = pickle.load(f)
			Np = data['Np']
			p_t = data['p_t']
			radius_dust = data['radius_dust']
			day = data['day']
			print(f"\nDays in the data: {day}\n{'-'*60}")
			for i in range(len(day)):
				data_shape = np.shape(p_t[i])
				print(f"\ntime = day **{day[i]}**,  Np = **{Np[i]}**")
				print(f"shape of data: {data_shape}")
				if data_shape[0] > 5:
					print(f"First 5 rows:\n{p_t[i][:5]}")
					print(f"Last 5 rows:\n{p_t[i][-5:]}")
				else: print(f"All rows:\n{p_t[i]}")
		except Exception as e:
			print(f"An error occurred while processing **{os.path.basename(datapath)}**: {e}")
	return


#------------------------------------------------------------------------------
if __name__ == "__main__":
	#quick_plot_particles()
	#check_pkl_structure("/home/linfel/linfel_scratch/rebound_exp/data_high_longterm_run_034-050/snapshot_data")
	check_metadata()
