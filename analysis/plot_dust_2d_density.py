import numpy as np
import os
import matplotlib.pyplot as plt
from ReadParticle import read_particle_frames

# data path
data_rootdir = "/u/lli22/ejecta_size_exp"
filenames = [os.path.join(data_rootdir, f"run_{i:03d}", "particles.txt") for i in range(1, 22)]

# Ensure output directory exists for saving frames
output_dir = "/u/lli22/ejecta_size_exp/particle_density"
os.makedirs(output_dir, exist_ok=True)

# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                          [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                          [-0.104839674791979,  0.124915784491013,  0.986612735258626]])

# sky north vector and convert it to 'Didymos System Barycenter' frame
sky_north = np.array([0, 0, 1])
r_sky_north = np.dot(SBC_rotate_DSB, sky_north.T)

# pixel dimension parameters
axlims = [-10000e3, 35000e3, -6000e3, 2000e3] # axis range [x_min, x_max, y_min, y_max] (m)
nx = 1125  # number of bins in x axis
ny = 200   # number of bins in y axis
xedges = np.linspace(axlims[0], axlims[1], nx + 1)
yedges = np.linspace(axlims[2], axlims[3], ny + 1)

# time slice parameters
target_day = 131.29
target_seconds = target_day * 86400.

for file in filenames:
    # Read particle data
    Np_seq, time, radii_dust, data_p = read_particle_frames(file)

    t_idx = np.searchsorted(time, target_seconds, side="left")

    p_t = data_p[t_idx]
    day = time[t_idx] / 86400.

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

    # Create a mask for points within the desired range
    mask = (
        (x_proj >= axlims[0]) & (x_proj <= axlims[1]) &
        (y_proj >= axlims[2]) & (y_proj <= axlims[3])
    )

    # Apply the mask
    x_filtered = x_proj[mask]
    y_filtered = y_proj[mask]

    # Create 2D histogram (density map)
    H, _, _ = np.histogram2d(
        x_filtered,  # x coordinate
        y_filtered,  # y coordinate
        bins=[xedges, yedges],
        density=False  # Set True if you want normalized density
    )

    # Plot the 2D histogram
    plt.figure(figsize=(8, 6))
    im = plt.imshow(
        H.T,
        origin='lower',
        extent=[val / 1e3 for val in axlims],
        aspect='equal',
        cmap='cividis'
    )
    plt.xlabel('Projected X [km]')
    plt.ylabel('Projected Y [km]')
    plt.title('2D Particle Density on View Plane')
 
    cbar = plt.colorbar(im, orientation='horizontal', pad=0.1)  # pad adjusts spacing
    cbar.set_label('Counts per bin')
		
    plt.grid(False)
    plt.tight_layout()

    output_name = os.path.join(output_dir, f"{radii_dust:.6f}.png")
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    plt.close()
