import numpy as np
import struct

Np_MAX_SAFE = 1000000
TIME_MAX_SAFE = 500 * 86400.

#--------------------------------------------------------------------------
def read_particle(file_particle, file_collide):
	# Initialize lists to hold data
	Np_seq = []  # time sequence of the number of particles 
	time = []
	data_p = []
	r_dust = None

	# Read particle file
	with open(file_particle, 'rb') as file:
		while True:
			try:
				Np_seq_new = struct.unpack('i', file.read(4))[0]  # Read int
				time_new = struct.unpack('d', file.read(8))[0]  # Read double
				r_dust = struct.unpack('d', file.read(8))[0]  # Read double
				data_new = np.fromfile(file, dtype=np.double, count=7*Np_seq_new).reshape((Np_seq_new, 7))
			except struct.error:
				break  # Break the loop if we run out of data to read

			Np_seq.append(Np_seq_new)
			time.append(time_new)
			data_p.append(data_new)
	
	Np_seq = np.array(Np_seq)
	time = np.array(time)
	Np_tot = Np_seq[0]

	# Process particle data into individual particle information at initial time step
	particle_all = [{'ID':         data_p[0][n, 0],
									 'id_collide': 0,
									 'time':       time[-1],
									 'pos':        [data_p[0][n, 1:4]],  # list tracking position, convert to np.array in the end
									 'vel':        [data_p[0][n, 4:7]]}  # list tracking velocity, convert to np.array in the end
									for n in range(Np_tot)]

	# Populate further time steps
	for i in range(1, len(time)):
		for n in range(Np_seq[i]):
			p_id = int(data_p[i][n, 0])
			particle_all[p_id-1]['pos'].append(data_p[i][n, 1:4])
			particle_all[p_id-1]['vel'].append(data_p[i][n, 4:7])
		print(f'No. {i+1} of output frame has been analyzed ({100.0 * (i+1) / len(Np_seq):.1f}%)', end="\r", flush=True)
	print("\noptimizing data structure ...", flush=True)
	for n in range(Np_tot):  # convert list to np.array
		particle_all[n]['pos'] = np.array(particle_all[n]['pos'])
		particle_all[n]['vel'] = np.array(particle_all[n]['vel'])

	# Read collision file
	data_c = []
	with open(file_collide, 'rb') as file:
		while True:
			try:
				pieces = struct.unpack('iid', file.read(16))  # Read int, int, double
			except struct.error:
				break  # Break the loop if we run out of data to read
			data_c.append(list(pieces))
    
	# Analyze collision data
	print("analyzing collision data ...", flush=True)
	N_colDidy, N_colDimor, N_escape = np.empty((0,2)), np.empty((0,2)), np.empty((0,2))
	for flag_remove, p_id, time_new in data_c:
		particle_all[p_id-1]['id_collide'] = flag_remove
		particle_all[p_id-1]['time'] = time_new  # time of removal
		if flag_remove == 1:
			N_colDidy = np.vstack(( N_colDidy, np.array([time_new, len(N_colDidy) + 1]) ))
			#N_colDidy.append((time_new, len(N_colDidy) + 1))
		elif flag_remove == 2:
			N_colDimor = np.vstack(( N_colDimor, np.array([time_new, len(N_colDimor) + 1]) ))
			#N_colDimor.append((time_new, len(N_colDimor) + 1))
		elif flag_remove == 3:
			N_escape = np.vstack(( N_escape, np.array([time_new, len(N_escape) + 1]) ))
			#N_escape.append((time_new, len(N_escape) + 1))
	
	print("processing completed!")
	return Np_seq, time, N_colDidy, N_colDimor, N_escape, r_dust, particle_all, data_c, data_p

#--------------------------------------------------------------------------
def read_particle_frames(file_particle):
	"""
	Reads all frames in the particle file
	"""
	# Initialize lists to hold data
	Np_seq = []  # time sequence of the number of particles 
	time = []
	data_p = []
	r_dust = None

	# Read particle file
	with open(file_particle, 'rb') as file:
		while True:
			try:
				Np_seq_new = struct.unpack('i', file.read(4))[0]  # Read int
				time_new = struct.unpack('d', file.read(8))[0]  # Read double
				r_dust = struct.unpack('d', file.read(8))[0]  # Read double
				data_new = np.fromfile(file, dtype=np.double, count=7*Np_seq_new).reshape((Np_seq_new, 7))
			except struct.error:
				break  # Break the loop if we run out of data to read

			Np_seq.append(Np_seq_new)
			time.append(time_new)
			data_p.append(data_new)
	
	Np_seq = np.array(Np_seq)
	time = np.array(time)
	Np_tot = Np_seq[0]

	print(f"{file_particle}  processing completed!", flush=True)
	return Np_seq, time, r_dust, data_p

#--------------------------------------------------------------------------
def read_specific_frame(file_particle, target_seconds):
	"""
	Reads a binary particle file and efficiently extracts the data
	for the first time frame that is at or after target_seconds.

	Args:
			file_particle (str): Path to the binary file.
			target_seconds (float or list[float]): One or more target times (seconds).

	Returns:
			dict: {target_time: (Np, time, r_dust, data)} for each target,
			or None if no suitable frame is found.
	"""
	target_seconds = sorted(target_seconds)  # ensure ascending order
	
	Np, time, data_p, r_dust = [], [], [], None
	tidx = 0
	n_targets = len(target_seconds)
	
	with open(file_particle, 'rb') as file:
		print(f"Processing particle data:  {file_particle}", flush=True)
		while tidx < n_targets:
			try:
				# Read the metadata for the current frame
				Np_current = struct.unpack('i', file.read(4))[0]
				time_current = struct.unpack('d', file.read(8))[0]
				if Np_current < 0 or Np_current > Np_MAX_SAFE:
					raise ValueError(f"Corrupt Np value ({Np_current}) detected.")
				if time_current < 0 or time_current > TIME_MAX_SAFE:
					raise ValueError(f"Corrupt time value ({time_current:.2f}) detected.")

				# If this frame's time is at or after our target...
				if time_current >= target_seconds[tidx]:
					# This is our frame. Read the rest of its data.
					r_dust = struct.unpack('d', file.read(8))[0]
					data_current = np.fromfile(file, dtype=np.double, count=7*Np_current).reshape((Np_current, 7))
					
					# Save result to lists and move to next target
					Np.append(Np_current)
					time.append(time_current)
					data_p.append(data_current)
					tidx += 1
					print(f"Time frame t={time_current} (day={time_current/86400:.2f}) processed.", flush=True)
				
				else:
					# Skip the large data block.
					bytes_to_skip = 8 + 7 * Np_current * 8
					file.seek(bytes_to_skip, 1) # 1 = seek from current position

			except (struct.error, EOFError, OSError, ValueError) as e:
				# This handles cases where the file ends unexpectedly.
				print(f"Warning: File terminated or seek failed due to truncation. Error: {e}", flush=True)
				break
	
	Np   = np.array(Np)
	time = np.array(time)
	return Np, time, r_dust, data_p

#--------------------------------------------------------------------------
def read_meta(file_particle):
	"""
	Read only metadata of the particle binary data.
	"""
	Np, time = [], []
	with open(file_particle, 'rb') as file:
		while True:
			try:
				Np_current = struct.unpack('i', file.read(4))[0]  # Read int
				time_current = struct.unpack('d', file.read(8))[0]  # Read double
				Np.append(Np_current)
				time.append(time_current)

				# skip the data block
				bytes_to_skip = 8 + 7 * Np_current * 8
				file.seek(bytes_to_skip, 1) # 1 = seek from current position

			except struct.error:
				break  # Break the loop if we run out of data to read
	
	Np   = np.array(Np)
	time = np.array(time)
	return Np, time
