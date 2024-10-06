import numpy as np
import matplotlib.pyplot as plt
from ReadParticle import read_particle 

# Constants
mass_system = 5.5e11
mu_system = 6.6743e-11 * 5.5e11
vol_didy = 0.2295409644951028e9
vol_dimor = 0.001830603200702610e9
mass_didy = mass_system * vol_didy / (vol_didy + vol_dimor)
mass_dimor = mass_system * vol_dimor / (vol_didy + vol_dimor)

# Read particle and collision data
Np_seq, time, N_colDidy, N_colDimor, N_escape, r_dust, particle, data_c, data_p = read_particle('particles.txt', 'collide.txt')

if (0):
    # Particles excluding escapers
    Np_seq_cor = Np_seq.copy()
    for i in range(len(N_escape)):
        index_t = np.where(np.array(time) >= N_escape[i,0],1,0)
        Np_seq_cor = Np_seq_cor + index_t    # uncollided particles 

    R_hill = 70000  # Hill radius, m
    time_hill = []  # Time reach the boundary of the hill radius, s
    N_hill = []  # Particle ID

    for i in range(3, len(particle)):
        for j in range(len(particle[i]['pos'])):
            dis_Didy = np.linalg.norm(particle[i]['pos'][j] - particle[0]['pos'][j])
            if dis_Didy > R_hill:
                time_hill.append(time[j])
                N_hill.append(particle[i]['ID'])
                if particle[i]['id_collide'] < 3:
                    data_c = np.vstack((data_c, [3, particle[i]['ID'], time[j]]))
                break

    time_sort = np.sort(time_hill)
    N_sort = np.arange(1, len(time_sort) + 1)

    # Particles including new-defined escapers
    Np_seq_cor_esc = Np_seq_cor.copy()
    for i in range(len(N_sort)):
        index_t = np.where(np.array(time) >= time_sort[i],1,0)
        Np_seq_cor_esc = Np_seq_cor_esc - index_t


# Plotting ejecta type
if (0):
    print("# Plotting ejecta type ...",flush=True)
    plt.figure()
    plt.plot(np.append(np.insert(N_colDidy[:,0],0,0),time[-1]) / 24 / 3600,
             np.append(np.insert(N_colDidy[:,1],0,0),N_colDidy[-1,1]) / Np_seq[0] * 100,
             label='Didymos collider', linewidth=1.5)
    plt.plot(np.append(np.insert(N_colDimor[:,0],0,0),time[-1]) / 24 / 3600,
             np.append(np.insert(N_colDimor[:,1],0,0),N_colDimor[-1,1]) / Np_seq[0] * 100,
             label='Dimorphos collider', linewidth=1.5)
    plt.plot(np.append(np.insert(N_escape[:,0],0,0),time[-1]) / 24 / 3600,
             np.append(np.insert(N_escape[:,1],0,0),N_escape[-1,1]) / Np_seq[0] * 100,
             label='Escaped ejecta', linewidth=1.5)
    #plt.plot(np.insert(time_sort, 0, 0) / 24 / 3600, np.insert(N_sort, 0, 0) / Np_seq[0] * 100, label='Escaped ejecta', linewidth=1.5)
    plt.plot(time / 24 / 3600,
             Np_seq / Np_seq[0] * 100,
             label='Remaining ejecta', linewidth=1.5)
    plt.xlabel('Time [days]')
    plt.ylabel('Ejecta type percentage [%]')
    plt.title('Dust particle radius r = 1 mm')
    plt.legend()
    plt.grid()
    plt.savefig('ejecta_type.png',dpi=300)
    plt.close()
    print("# Ejecta type completed!\n",flush=True)

# Dimorphos orbit
if (0):
    print("# Plotting Dimorphos orbit ...",flush=True)
    didy_pos  = particle[0]['pos']
    dimor_pos = particle[1]['pos']
    n_p = 10
    np_pos    = particle[n_p]['pos']

    plt.figure().add_subplot(projection='3d')
    plt.plot(dimor_pos[:, 0] - didy_pos[:, 0],
             dimor_pos[:, 1] - didy_pos[:, 1],
             dimor_pos[:, 2] - didy_pos[:, 2], 
             label = 'Dimorphos orbit', linewidth=1.5, color='k')
    #plt.plot(np_pos[:, 0] - didy_pos[0:len(np_pos), 0],
    #         np_pos[:, 1] - didy_pos[0:len(np_pos), 1],
    #         np_pos[:, 2] - didy_pos[0:len(np_pos), 2],
    #         label = f'particle {n_p} orbit', linewidth=1.5)
    plt.savefig('dimor_orbit.png',dpi=300)
    plt.close()
    print("# Dimorphos orbit completed!\n",flush=True)


# Scatter plot for time step N_t
if (0):
    print("# Plotting locations of particles at a time slice ...",flush=True)
    N_t = 0
    plt.figure()
    # plot Didymos and Dimorphos
    plt.scatter(data_p[N_t][0, 1], data_p[N_t][0, 2], c='red', s=10)
    plt.scatter(data_p[N_t][1, 1], data_p[N_t][1, 2], c='red', s=8)
    # plot dust particles
    dust_sc = plt.scatter(data_p[N_t][3:, 1], data_p[N_t][3:, 2], c=data_p[N_t][3:, 3], s=2, cmap='viridis')
    cb = plt.colorbar(dust_sc)
    cb.set_label(f'z [m]', fontsize=12)
    
    #plt.axes().set_aspect('equal')
    plt.title(f't = {time[N_t]:.1f} s   Dust particle radius r = 1 mm')
    plt.grid()
    plt.savefig(f't{N_t}_scatter.png',dpi=300)
    plt.close()
    print("# Locations of particles completed!\n",flush=True)

# read particle a_e data
p_ae = np.genfromtxt("a_e.csv", delimiter=',', skip_header=1)

# Scatter plot for dust fate
if (0):
    print("# Plotting dust fate scatter plot ...",flush=True)
    fate_list = np.array([p['id_collide'] for p in particle[3:]])
    plt.figure()
    plt.scatter(p_ae[fate_list == 0,2], p_ae[fate_list == 0,1], s=5, label='ID = 0')
    plt.scatter(p_ae[fate_list == 1,2], p_ae[fate_list == 1,1], s=5, label='ID = 1')
    plt.scatter(p_ae[fate_list == 2,2], p_ae[fate_list == 2,1], s=5, label='ID = 2')
    plt.scatter(p_ae[fate_list == 3,2], p_ae[fate_list == 3,1], s=5, label='ID = 3')
    plt.title('Dust particle radius r = 1 mm')
    plt.xlabel('Eccentricity')
    plt.ylabel('Semimajor axis [m]')
    plt.legend()
    plt.savefig('a_e.png',dpi=300)
    plt.close()
    print("# Dust fate scatter plot completed!\n",flush=True)

# get IDs based on particle fate classification
remain_id, didy_col_id, dimor_col_id, esc_id = [], [], [], []
for i in range(len(particle)):
    fate_flag = particle[i]['id_collide']
    ID = int(particle[i]['ID'])
    if fate_flag == 0:
        remain_id.append(ID)
    elif fate_flag == 1:
        didy_col_id.append(ID)
    elif fate_flag == 2:
        dimor_col_id.append(ID)
    elif fate_flag == 3:
        esc_id.append(ID)
    else:
        print('collision category error',flush=True)
    print(f'fate of No. {i+1} particle classified. ({100.0 * (i+1) / len(particle):.1f}%)', end="\r", flush=True)
print("\n", flush=True)

remain_id = np.array(remain_id)
didy_col_id = np.array(didy_col_id)
dimor_col_id = np.array(dimor_col_id)
esc_id = np.array(esc_id)

# Histogram for semi-major axis and eccentricity
## semimajor axis
if (0):
    print("# Plotting histogram for semimajor axis ...",flush=True)
    fig,axs = plt.subplots(2,2, figsize=(8,6))
    bins = np.linspace(-100,11000,112)
    
    if len(didy_col_id) > 0:
        counts, bin_edges = np.histogram(p_ae[didy_col_id - 4,1], bins=bins)
        axs[0,0].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Didymos collider')
        axs[0,0].set_yscale('log')
        axs[0,0].legend(fontsize='small')
    
    if len(dimor_col_id) > 0:
        counts, bin_edges = np.histogram(p_ae[dimor_col_id - 4,1], bins=bins)
        axs[0,1].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Dimorphos collider')
        axs[0,1].set_yscale('log')
        axs[0,1].legend(fontsize='small')
    
    if len(esc_id) > 0:
        counts, bin_edges = np.histogram(p_ae[esc_id - 4,1], bins=bins)
        axs[1,0].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Escaped ejecta')
        axs[1,0].set_yscale('log')
        axs[1,0].legend(fontsize='small')
    
    if len(remain_id) > 3: # exclude Didymos, Dimorphos and Sun //need improvement to capture Didy Dimor Sun
        counts, bin_edges = np.histogram(p_ae[remain_id - 4,1], bins=bins)
        axs[1,1].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Remaining ejecta')
        axs[1,1].set_yscale('log')
        axs[1,1].legend(fontsize='small')

    fig.text(0.5, 0.04, 'Semimajor axis [m]' ,ha='center', va='center')
    fig.text(0.04, 0.5, 'Number', ha='center', va='center', rotation='vertical')
    plt.suptitle('Dust particle radius r = 1 mm')
    fig.savefig('semimajor_axis.png',dpi=300)
    plt.close()
    print("# Histogram for semimajor axis completed!\n",flush=True)

## eccentricity
if (0):
    print("# Plotting histogram for eccentricity ...",flush=True)
    fig,axs = plt.subplots(2,2, figsize=(8,6))
    bins = np.linspace(-1,5,150)
    
    if len(didy_col_id) > 0:
        counts, bin_edges = np.histogram(p_ae[didy_col_id - 4,2], bins=bins)
        axs[0,0].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Didymos collider')
        axs[0,0].set_yscale('log')
        axs[0,0].legend(fontsize='small')
    
    if len(dimor_col_id) > 0:
        counts, bin_edges = np.histogram(p_ae[dimor_col_id - 4,2], bins=bins)
        axs[0,1].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Dimorphos collider')
        axs[0,1].set_yscale('log')
        axs[0,1].legend(fontsize='small')
    
    if len(esc_id) > 0:
        counts, bin_edges = np.histogram(p_ae[esc_id - 4,2], bins=bins)
        axs[1,0].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Escaped ejecta')
        axs[1,0].set_yscale('log')
        axs[1,0].legend(fontsize='small')
    
    if len(remain_id) > 3: # exclude Didymos, Dimorphos and Sun //need improvement to capture Didy Dimor Sun
        counts, bin_edges = np.histogram(p_ae[remain_id - 4,2], bins=bins)
        axs[1,1].bar(bin_edges[:-1], counts, width=bin_edges[1]-bin_edges[0], align='edge', label='Remaining ejecta')
        axs[1,1].set_yscale('log')
        axs[1,1].legend(fontsize='small')

    fig.text(0.5, 0.04, 'Eccentricity' ,ha='center', va='center')
    fig.text(0.04, 0.5, 'Number', ha='center', va='center', rotation='vertical')
    plt.suptitle('Dust particle radius r = 1 mm')
    fig.savefig('eccentricity.png',dpi=300)
    plt.close()
    print("# Histogram for eccentricity completed!\n",flush=True)
