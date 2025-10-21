import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from ReadParticle import read_particle_frames
import miepython as mie

# data path
data_rootdir = "/home/linfel/linfel_turbo/rebound_exp/snapshot_data/day11.88"
filenames = [os.path.join(data_rootdir, f"particle_{i:03d}_day11.88.pkl") for i in [1,11,21,31,40,50]]

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/linfel_turbo/rebound_exp/plots/plots_day11.88"
os.makedirs(output_dir, exist_ok=True)

# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                          [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                          [-0.104839674791979,  0.124915784491013,  0.986612735258626]])

# calculate sky north vector in 'Sun Body Center' and 
# convert it to 'Didymos System Barycenter' frame
obliq_earth = np.deg2rad(23.4392911)
sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
r_sky_north = SBC_rotate_DSB @ sky_north

# optical parameters
m = 1.5           # refractive index of particle
lambda0 = 500e-9  # wavelength in vacuum (m)

# pixel dimension parameters
#axlims = [-28000e3, 28000e3, -28000e3, 28000e3] # axis range [x_min, x_max, y_min, y_max] (m)
axlims = [-600e3, 600e3, -600e3, 600e3]
nx = 4000  # number of bins in x axis
ny = 4000   # number of bins in y axis
xedges = np.linspace(axlims[0], axlims[1], nx + 1)
yedges = np.linspace(axlims[2], axlims[3], ny + 1)

# initialize total intensity as 1e-20
total_inten = np.zeros((nx, ny)) + 1e-20

for file in filenames:
    # Read particle data
    print(f"processing {file} ......", flush=True)
    with open(file, 'rb') as f:
        data = pickle.load(f)
        p_t = data['p_t']
        radius_dust = data['radius_dust']
        day = data['day']

    # position vector
    r_sun = p_t[2, 1:4]    # [x, y, z] of the Sun
    r_earth = p_t[3, 1:4]  # [x, y, z] of the Earth   (as a proxy of HST)
    r_dust = p_t[4:, 1:4]  # [x, y, z] of all dust particles

    # calculate the two basis vectors of the projection plane
    l1 = np.cross(r_sky_north, r_earth)
    l2 = np.cross(r_earth, l1)
    l1 = l1 / np.linalg.norm(l1)
    l2 = l2 / np.linalg.norm(l2)

    # create an array to record coordinates of particles projected onto the plane
    p_projected = np.zeros((len(p_t),3), dtype=float)
    p_projected[:, 0] = p_t[:, 0]                  # copy the column of particle ID
    p_projected[:, 1] = np.dot(p_t[:, 1:4], l1)    # x coordinate
    p_projected[:, 2] = np.dot(p_t[:, 1:4], l2)    # y coordinate

    """2D density histogram"""
    hist_data = p_projected[4:]

    # Extract x and y projected coordinates
    x_proj = hist_data[:, 1]
    y_proj = hist_data[:, 2]

    # Mask out particles outside of range
    mask = (
        (x_proj >= axlims[0]) & (x_proj <= axlims[1]) &
        (y_proj >= axlims[2]) & (y_proj <= axlims[3])
    )
    x_filtered = x_proj[mask]
    y_filtered = y_proj[mask]

    # Create 2D histogram
    px_den, _, _ = np.histogram2d( # particle density of each pixel
        x_filtered,        # x coordinate
        y_filtered,        # y coordinate
        bins=[xedges, yedges],
        density=False  # Set True if you want normalized density
    )
    px_den = px_den.T

    # set weight
    W = 1e6
    #weight = 10 ** ((np.log10(radius_dust)+1) * np.log10(W) / (-3))    
    weight = 1

    # scattering intensity (assume scattering phase angle constant for all particles)
    qext, qsca, qback, g = mie.efficiencies(m, 2*radius_dust, lambda0)
    p_func = 1
    px_inten = qsca * (np.pi * radius_dust**2) * px_den * weight * p_func

    # Accumulate intensity of each pixel
    total_inten += px_inten

# take log10 of intensity
log10_inten = np.log10(total_inten)

# Create meshgrid for bin edges
X, Y = np.meshgrid(xedges/1e3, yedges/1e3)   # km

# Plot using pcolor
plt.figure(figsize=(10, 8))
pc = plt.pcolormesh(X, Y, log10_inten, cmap='cividis', shading='auto', vmin=-10, vmax=np.nanmax(log10_inten))

plt.xlabel('Projected X [km]')
plt.ylabel('Projected Y [km]')
plt.title('Intensity on View Plane (Custom Weights)')
plt.grid(True, linestyle='--', linewidth=0.1, color='red', alpha=0.7)
plt.gca().set_aspect('equal', adjustable='box')

cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.5, aspect=30)  # pad adjusts spacing
cbar.set_label(r'$\log_{10}$(Nondimensional Intensity)')

plt.tight_layout()

output_name = os.path.join(output_dir, f"day11.88_intensity_1.png")
#output_name = os.path.join(output_dir, f"tail_synt_wt1.png")
plt.savefig(output_name, dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()
