import numpy as np
import struct

def read_particle(file_particle, file_collide):
    # Initialize lists to hold data
    N_particle = []
    time = []
    data_p = []
    r_dust = None

    # Read particle file
    with open(file_particle, 'rb') as file:
        while True:
            try:
                N_particle_new = struct.unpack('i', file.read(4))[0]  # Read int
                time_new = struct.unpack('d', file.read(8))[0]  # Read double
                r_dust = struct.unpack('d', file.read(8))[0]  # Read double
                data_new = np.fromfile(file, dtype=np.double, count=7*N_particle_new).reshape((N_particle_new, 7))
            except struct.error:
                break  # Break the loop if we run out of data to read

            N_particle.append(N_particle_new)
            time.append(time_new)
            data_p.append(data_new)

    # Read collision file
    data_c = []
    with open(file_collide, 'rb') as file:
        while True:
            try:
                pieces = struct.unpack('iid', file.read(16))  # Read int, int, double
            except struct.error:
                break  # Break the loop if we run out of data to read
            data_c.append(pieces)

    # Process data into individual particle information
    particle_all = [{'ID': data_p[0][n, 0], 'id_collide': 0, 'time': time[-1], 'pos': data_p[0][n, 1:4], 'vel': data_p[0][n, 4:7]}
                    for n in range(N_particle[0])]

    # Populate further time steps
    for i in range(1, len(N_particle)):
        print(f'No. {i} of output frame has been analyzed ({100.0 * i / len(N_particle):.1f}%)')
        for n in range(N_particle[i]):
            p_id = int(data_p[i][n, 0])
            particle_all[p_id]['pos'] = np.vstack((particle_all[p_id]['pos'], data_p[i][n, 1:4]))
            particle_all[p_id]['vel'] = np.vstack((particle_all[p_id]['vel'], data_p[i][n, 4:7]))

    # Analyze collision data
    N_colDidy, N_colDimor, N_escape = [], [], []
    for flag_remove, p_id, time_new in data_c:
        particle_all[p_id]['id_collide'] = flag_remove
        particle_all[p_id]['time'] = time_new
        if flag_remove == 1:
            N_colDidy.append((time_new, len(N_colDidy) + 1))
        elif flag_remove == 2:
            N_colDimor.append((time_new, len(N_colDimor) + 1))
        elif flag_remove == 3:
            N_escape.append((time_new, len(N_escape) + 1))

    return N_particle, time, N_colDidy, N_colDimor, N_escape, r_dust, particle_all, data_c, data_p

