import os
from ReadParticle import read_particle_frames
import numpy as np

DATA_ROOTDIR = "/home/linfel/linfel_turbo/rebound_exp/data_high_shortterm_run"
RUN_NUMBERS = range(34,51)

MIN_ALLOWED_RADIUS = 0.0
MAX_ALLOWED_VALUE = 2 * 1.496e+11  # 2 AU
MAX_ALLOWED_NUMBER = 1000000       # max allowed number of particles

for run_idx in RUN_NUMBERS:
	print(f"\nprocessing run_{run_idx:03d}")
	data_path = os.path.join(DATA_ROOTDIR, f"run_{run_idx:03d}", "particles.txt")
	Np, time, r_dust, data_p = read_particle_frames(data_path)

	print(f"r={r_dust}")
	if r_dust < MIN_ALLOWED_RADIUS:
		print(f"**Bad Data Warning:** Negative dust radius found in run_{run_idx:03d}!")

	for i in range(len(time)):
		current_day = time[i] / 86400.
		current_data = data_p[i]
		data_shape = np.shape(current_data)
		
		print(f"\nday: {current_day},    Np={Np[i]}")
		print(f"shape of data: {data_shape}")

		# --- Data Integrity Checks ---
		assert data_shape[0] == Np[i], "**Bad Data Alert:** Number of rows doesn't match number of particles"
		assert Np[i] < MAX_ALLOWED_NUMBER, "**Bad Data Alert:** Invalid number of particles found at day {current_day} in run_{run_idx:03d}!"

		# 1. Check for NaN or Inf (Not a Number or Infinity)
		if np.any(np.isnan(current_data)):
			print(f"**Bad Data Alert:** NaN value found at day {current_day} in run_{run_idx:03d}!")
		
		if np.any(np.isinf(current_data)):
			print(f"**Bad Data Alert:** Infinity (Inf) value found at day {current_day} in run_{run_idx:03d}!")

		# 2. Check for Unphysically Large Magnitudes (Explosion/Escape)
		# Check if the absolute value of any coordinate/velocity exceeds a threshold
		if np.any(np.abs(current_data) > MAX_ALLOWED_VALUE):
			print(f"**Bad Data Warning:** Extreme position/velocity value (> {MAX_ALLOWED_VALUE:.2E}) found at day {current_day} in run_{run_idx:03d}!")
			print(f"Max absolute value found: {np.max(np.abs(current_data))}")
