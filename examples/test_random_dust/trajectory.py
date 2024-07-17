import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import norm
from ReadParticle import read_particle 

# Constants
mass_system = 5.5e11
mu_system = 6.6743e-11 * 5.5e11
vol_didy = 0.2295409644951028e9
vol_dimor = 0.001830603200702610e9
mass_didy = mass_system * vol_didy / (vol_didy + vol_dimor)
mass_dimor = mass_system * vol_dimor / (vol_didy + vol_dimor)

# Read particle and collision data
N_particle, time, N_colDidy, N_colDimor, N_escape, r_dust, particle, data_c, data_p = read_particle('particles.txt', 'collide.txt')

# Particles excluding escapers
N_particle_cor = N_particle.copy()
for i in range(len(N_escape)):
    index_t = np.where(time >= N_escape[i, 0])[0]
    if len(index_t) > 0:
        N_particle_cor[index_t[0]:] += 1

R_hill = 70000  # Hill radius, m
time_hill = []  # Time reach the boundary of the hill radius, s
N_hill = []  # Particle ID

ii = 1
for i in range(4, len(particle)):
    for j in range(len(particle[i]['pos'])):
        dis_Didy = np.linalg.norm(particle[i]['pos'][j] - particle[1]['pos'][j])
        if dis_Didy > R_hill:
            time_hill.append(time[j])
            N_hill.append(particle[i]['ID'])
            ii += 1
            if particle[i]['id_collide'] < 3:
                data_c = np.vstack((data_c, [3, particle[i]['ID'], time[j]]))
            break

time_sort = np.sort(time_hill)
N_sort = np.arange(1, len(time_sort) + 1)

# Particles including new-defined escapers
N_particle_cor_esc = N_particle_cor.copy()
for i in range(len(N_sort)):
    index_t = np.where(time >= time_sort[i])[0]
    if len(index_t) > 0:
        N_particle_cor_esc[index_t[0]:] -= 1

# Plotting
plt.figure()
plt.plot(np.insert(N_colDidy[:, 0], 0, 0) / 24 / 3600, np.insert(N_colDidy[:, 1], 0, 0) / N_particle[0] * 100, label='Didymos collider', linewidth=1.5)
plt.plot(np.insert(N_colDimor[:, 0], 0, 0) / 24 / 3600, np.insert(N_colDimor[:, 1], 0, 0) / N_particle[0] * 100, label='Dimorphos collider', linewidth=1.5)
plt.plot(np.insert(time_sort, 0, 0) / 24 / 3600, np.insert(N_sort, 0, 0) / N_particle[0] * 100, label='Escaped ejecta', linewidth=1.5)
plt.plot(time / 24 / 3600, N_particle_cor_esc / N_particle[0] * 100, label='Remaining ejecta', linewidth=1.5)
plt.xlabel('Time [days]')
plt.ylabel('Ejecta type percentage [%]')
plt.title('Dust particle radius r = 1 mm')
plt.legend()
plt.grid()
plt.savefig('ejecta_type.png',dpi=300)
plt.close()

# Dimorphos orbit
plt.figure()
plt.plot(particle[2]['pos'][:, 0] - particle[1]['pos'][:, 0],
         particle[2]['pos'][:, 1] - particle[1]['pos'][:, 1],
         particle[2]['pos'][:, 2] - particle[1]['pos'][:, 2], linewidth=2.0, color='k')
N_p = 200
plt.plot(particle[N_p]['pos'][:, 0] - particle[1]['pos'][0:len(particle[N_p]['pos']), 0],
         particle[N_p]['pos'][:, 1] - particle[1]['pos'][0:len(particle[N_p]['pos']), 1],
         particle[N_p]['pos'][:, 2] - particle[1]['pos'][0:len(particle[N_p]['pos']), 2])
plt.savefig('dimor_orbit.png',dpi=300)
plt.close()

# Repeat previous plot for consistency
#plt.figure()
#plt.plot(np.insert(N_colDidy[:, 0], 0, 0) / 24 / 3600, np.insert(N_colDidy[:, 1], 0, 0) / N_particle[0] * 100, label='Didymos collider', linewidth=1.5)
#plt.plot(np.insert(N_colDimor[:, 0], 0, 0) / 24 / 3600, np.insert(N_colDimor[:, 1], 0, 0) / N_particle[0] * 100, label='Dimorphos collider', linewidth=1.5)
#plt.plot(np.insert(N_escape[:, 0], 0, 0) / 24 / 3600, np.insert(N_escape[:, 1], 0, 0) / N_particle[0] * 100, label='Escaped ejecta', linewidth=1.5)
#plt.plot(time / 24 / 3600, N_particle / N_particle[0] * 100, label='Remaining ejecta', linewidth=1.5)
#plt.xlabel('Time [days]')
#plt.ylabel('Ejecta type percentage [%]')
#plt.title('Dust particle radius r = 1 mm')
#plt.legend()
#plt.grid()
#plt.show()

# Histogram for semi-major axis and eccentricity
a_p = np.zeros(300*20)
e_p = np.zeros(300*20)
for i in range(1, 301):
    for j in range(20):
        a_p[j + (i - 1) * 20] = 500.0 + (i - 1) * 20.0
        e_p[j + (i - 1) * 20] = (j - 1) / 20

plt.figure()
plt.hist(a_p[data_c[data_c[:, 0] == 1, 1] - 3], bins=np.arange(0, 6000, 200), label='Didymos collider')
plt.hist(a_p[data_c[data_c[:, 0] == 2, 1] - 3], bins=np.arange(0, 6000, 200), label='Dimorphos collider')
plt.hist(a_p[data_c[data_c[:, 0] == 3, 1] - 3], bins=np.arange(0, 6000, 200), label='Escaped ejecta')
plt.hist(a_p[data_p[-1][4:, 0] - 3], bins=np.arange(0, 6000, 200), label='Remaining ejecta')
plt.xlabel('Semimajor axis [m]')
plt.ylabel('Number')
plt.title('Dust particle radius r = 1 mm')
plt.legend()
plt.grid()
plt.savefig('semimajor_axis.png',dpi=300)
plt.close()

plt.figure()
plt.hist(e_p[data_c[data_c[:, 0] == 1, 1] - 3], bins=np.arange(0, 1, 0.05), label='Didymos collider')
plt.hist(e_p[data_c[data_c[:, 0] == 2, 1] - 3], bins=np.arange(0, 1, 0.05), label='Dimorphos collider')
plt.hist(e_p[data_c[data_c[:, 0] == 3, 1] - 3], bins=np.arange(0, 1, 0.05), label='Escaped ejecta')
plt.hist(e_p[data_p[-1][4:, 0] - 3], bins=np.arange(0, 1, 0.05), label='Remaining ejecta')
plt.xlabel('Eccentricity')
plt.ylabel('Number')
plt.title('Dust particle radius r = 1 mm')
plt.legend()
plt.grid()
plt.savefig('eccentricity.png',dpi=300)
plt.close()

# Scatter plot
color_code = np.array([p['id_collide'] for p in particle[3:]])
plt.figure()
plt.scatter(e_p[color_code == 0], a_p[color_code == 0], marker='o', linewidth=2.0, label='ID = 0')
plt.scatter(e_p[color_code == 1], a_p[color_code == 1], marker='v', linewidth=2.0, label='ID = 1')
plt.scatter(e_p[color_code == 2], a_p[color_code == 2], marker='p', linewidth=2.0, label='ID = 2')
plt.scatter(e_p[color_code == 3], a_p[color_code == 3], marker='d', linewidth=2.0, label='ID = 3')
plt.title('Dust particle radius r = 1 mm')
plt.legend()
plt.grid()
plt.savefig('a_e.png',dpi=300)
plt.close()

# Scatter plot for time step N_t
N_t = 12
plt.figure()
plt.scatter(data_p[N_t][3:, 2], data_p[N_t][3:, 3], c=data_p[N_t][3:, 1], s=10, cmap='viridis', edgecolor='k')
plt.colorbar()
plt.title('Dust particle radius r = 1 mm')
plt.grid()
plt.savefig('timeslice_scatter.png',dpi=300)
plt.close()

