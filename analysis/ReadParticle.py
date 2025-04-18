import numpy as np
import struct

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

def read_particle_frames(file_particle):
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
