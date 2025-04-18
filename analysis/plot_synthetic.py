import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from ReadParticle import read_particle_frames
import miepython as mie

# data path
data_rootdir = "/home/linfel/linfel_data/ejecta_size_exp/day131.29"
filenames = [os.path.join(data_rootdir, f"particle_{i:03d}_day131.29.pkl") for i in range(1, 3)]

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/linfel_data/ejecta_size_exp/synthetic_4_18"
os.makedirs(output_dir, exist_ok=True)

# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                          [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                          [-0.104839674791979,  0.124915784491013,  0.986612735258626]])

# sky north vector and convert it to 'Didymos System Barycenter' frame
sky_north = np.array([0, 0, 1])
r_sky_north = np.dot(SBC_rotate_DSB, sky_north.T)

# optical parameters
m = 1.5           # refractive index of particle
lambda0 = 500e-9  # wavelength in vacuum (m)

# pixel dimension parameters
axlims = [-10000e3, 35000e3, -6000e3, 2000e3] # axis range [x_min, x_max, y_min, y_max] (m)
nx = 1125  # number of bins in x axis
ny = 200   # number of bins in y axis
xedges = np.linspace(axlims[0], axlims[1], nx + 1)
yedges = np.linspace(axlims[2], axlims[3], ny + 1)

# initialize total intensity as 0
total_inten = np.zeros((nx, ny))

for file in filenames:
    # Read particle data
    with open(file, 'rb') as f:
        data = pickle.load(f)
        p_t = data['p_t']
        radii_dust = data['radii_dust']
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

    # set weight
    #W = 1e6
    #weight = 10 ** ((np.log10(radii_dust)+1) * np.log10(W) / (-3))    
    weight = 1

    # scattering intensity (assume scattering phase angle constant for all particles)
    qext, qsca, qback, g = mie.efficiencies(m, 2*radii_dust, lambda0)
    p_func = 1
    px_inten = qsca * (np.pi * radii_dust**2) * px_den * weight * p_func

    # Accumulate intensity of each pixel
    total_inten += px_inten

# Plot the 2D histogram
plt.figure(figsize=(12, 9))
im = plt.imshow(
    np.log(total_inten).T,
    origin='lower',
    extent=[val / 1e3 for val in axlims],
    aspect='equal',
    cmap='cividis'
)
plt.xlabel('Projected X [km]')
plt.ylabel('Projected Y [km]')
plt.title('2D Particle Density on View Plane')
plt.grid(False)

cbar = plt.colorbar(im, orientation='horizontal', pad=0.1)  # pad adjusts spacing
cbar.set_label('Nondimensional Intensity')

plt.tight_layout()

output_name = os.path.join(output_dir, f"tail_synt.png")
plt.savefig(output_name, dpi=300, bbox_inches='tight')
#plt.show()
plt.close()
