"""
Quick plot top view, HST view and side view
"""

import os
import numpy as np
from ReadParticle import read_particle_frames
from frame_renderer import init_worker, render_frame_topview, render_frame_HSTview, render_frame_sideview

def quickview_3(data_path, output_dir = '.', target_days=None, axis_lim=3e3):
	# Read particle data
	Np_seq, time, r_dust, data_p = read_particle_frames(data_path)

	# Ensure output directory exists for saving frames
	os.makedirs(output_dir, exist_ok=True)
	
	if target_days is None:
		target_seconds = time
	else:
		target_seconds = np.array(target_days) * 86400.0
    
	t_idx = np.searchsorted(time, target_seconds, side="left")
    
	if np.isscalar(axis_lim):
		axis_lim = [axis_lim] * len(t_idx)
        
	for i in range(len(t_idx)):
		init_worker(data_p, time, axis_lim[i], output_dir)

		render_frame_topview(t_idx[i])
		render_frame_sideview(t_idx[i])
		render_frame_HSTview(t_idx[i])
    
if __name__ == "__main__":
	DATA_ROOTDIR = "/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_run_BS"
	RUN_NUMBERS = range(1, 34) # Define the run numbers you want to process
	PLOT_DAYS = [0, 1.14, 5.7, 8.8, 11.86, 14.9]
	AXIS_LIM = [3e3, 200e3, 200e3, 200e3, 300e3, 300e3]        # correspond to each day

	for run_idx in RUN_NUMBERS:
		print(f"processing run_{run_idx:03d}")
		data_path = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")
		output_dir = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "postprocess/quickview")
		quickview_3(data_path, output_dir, PLOT_DAYS, AXIS_LIM)
