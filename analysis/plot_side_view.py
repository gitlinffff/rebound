import numpy as np
import os, pickle
import matplotlib.pyplot as plt
from ReadParticle import read_particle_frames

# data path
data_rootdir = "/home/linfel/rebound_exp"
filename = os.path.join(data_rootdir, f"particles.txt")

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/rebound_exp/output_singleframe"
os.makedirs(output_dir, exist_ok=True)

# Read particle data
Np_seq, time, radii_dust, data_p = read_particle_frames(filename)
p_t = data_p[0]
day = time[0]/86400

# position vector
r_dimor = p_t[1, 1:4]  # [x, y, z] of Dimorphos
v_dimor = p_t[1, 4:7]  # [vx, vy, vz] of Dimorphos

# calculate the two basis vectors of the projection plane
l2 = np.cross(v_dimor, r_dimor)
l1 = r_dimor / np.linalg.norm(r_dimor)
l2 = l2 / np.linalg.norm(l2)

# create an array to record coordinates of particles projected onto the plane
p_projected = np.zeros((len(p_t),3), dtype=float)
p_projected[:, 0] = p_t[:, 0]                  # copy the column of particle ID
p_projected[:, 1] = np.dot(p_t[:, 1:4], l1)    # x coordinate
p_projected[:, 2] = np.dot(p_t[:, 1:4], l2)    # y coordinate

# units to km
p_projected[:, 1:] = p_projected[:, 1:] / 1e3

# Create a new figure for this frame
fig, ax = plt.subplots(figsize=(8, 8))

# Determine axis limit
axlims = 3.  # km
if axlims is None:
		print("No axis limits specified, using default limits.")
elif isinstance(axlims, list):
		ax.set_xlim(axlims[0], axlims[1])
		ax.set_ylim(axlims[2], axlims[3])
else:
		ax.set_xlim(-axlims, axlims)
		ax.set_ylim(-axlims, axlims)

# plot Didymos and Dimorphos
ax.scatter(p_projected[0, 1], p_projected[0, 2], c='red', s=10, zorder=3, label='Didymos')
ax.scatter(p_projected[1, 1], p_projected[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

# Plot Sun direction relative to Didymos System Barycenter
sun_x, sun_y = p_projected[2, 1], p_projected[2, 2]
sun_distance = (sun_x**2 + sun_y**2) ** 0.5
arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')

# Plot Earth direction
earth_x, earth_y = p_projected[3, 1], p_projected[3, 2]
earth_distance = (earth_x**2 + earth_y**2) ** 0.5
arrow_x = earth_x / earth_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
arrow_y = earth_y / earth_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='green', label='Earth Direction')

ax.set_xlabel('X / km')
ax.set_ylabel('Z / km')
ax.grid()
ax.legend(loc='upper right')

# Set title
ax.set_title(f't = {day} days   Side view positions')

# Save the frame as a PNG image
plt.savefig(os.path.join(output_dir, f"sideview.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
