import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
import numpy as np

# Define global variables for shared data
shared_data = {}

def init_worker(data_p, time, axis_lim_list, output_dir):
    """
    Initialize global shared data for each worker process.
    This function is called once per worker during pool initialization.
    
    Args:
        data_p (list): A list of ndarray containing particle data for all frames.
        time (ndarray): Time data for all frames.
        axis_lim_list (list or float): List of axis limits for each frame (if rendering video) 
                                       or a single float value for one frame.
    """
    global shared_data
    shared_data['data_p'] = data_p
    shared_data['time'] = time
    shared_data['axis_lim_list'] = axis_lim_list
    shared_data['output_dir'] = output_dir

# Function to convert meters to kilometers for axis labels
def m_to_km(x, _):
    return f'{x / 1e3:.0f}'

def render_frame_topview(frame):
    """
    Renders a single frame of top view and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.
        output_dir (str): Directory to save the output image.

    Returns:
        str: The path of the saved frame image.
    """
    
    global shared_data
    data_p = shared_data['data_p']
    time = shared_data['time']
    axis_lim_list = shared_data['axis_lim_list']
    output_dir = shared_data['output_dir']

    print(f"Rendering frame {frame}...", end="\r", flush=True)
    p_t = data_p[frame]
    sec = time[frame]
    
    # Determine axis limit
    if isinstance(axis_lim_list, list):
        axis_lim = axis_lim_list[frame]  # Use frame-specific axis limit
    else:
        axis_lim = axis_lim_list  # Use single float value for a single frame
    
    # Create a new figure for this frame
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot Didymos and Dimorphos
    ax.scatter(p_t[0, 1], p_t[0, 2], c='red', s=10, zorder=3, label='Didymos')
    ax.scatter(p_t[1, 1], p_t[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

    # Plot dust particles
    speed = (p_t[4:, 4]**2 + p_t[4:, 5]**2 + p_t[4:, 6]**2)**0.5
    dust_sc = ax.scatter(p_t[4:, 1], p_t[4:, 2], c=speed, s=2, cmap='viridis')
    cb = fig.colorbar(dust_sc, ax=ax, shrink=0.8, aspect=20)
    cb.set_label(f'speed (m/s)', fontsize=12)

    # Plot Sun direction
    sun_x, sun_y = p_t[2, 1], p_t[2, 2]
    sun_distance = (sun_x**2 + sun_y**2) ** 0.5
    arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')

    # Plot Earth direction
    earth_x, earth_y = p_t[3, 1], p_t[3, 2]
    earth_distance = (earth_x**2 + earth_y**2) ** 0.5
    arrow_x = earth_x / earth_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    arrow_y = earth_y / earth_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='green', label='Earth Direction')

    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_xlabel('x / km')
    ax.set_ylabel('y / km')
    ax.grid()
    ax.legend()

    # Set title
    ax.set_title(f't = {sec/86400:.2f} days   Dust radius r = 1 mm')

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return frame_filename


def render_frame_HSTview(frame):
    """
    Renders a single frame of HST view and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.
        data_p (ndarray): Particle data for all frames.
        time (ndarray): Time data for all frames.
        axis_lim_list (list or float): List of axis limits for each frame (if rendering video) 
                                       or a single float value for one frame.
        output_dir (str): Directory to save the output image.

    Returns:
        str: The path of the saved frame image.
    """
    global shared_data
    data_p = shared_data['data_p']
    time = shared_data['time']
    axis_lim_list = shared_data['axis_lim_list']
    output_dir = shared_data['output_dir']
    
    print(f"Rendering frame {frame}...", end="\r", flush=True)
    p_t = data_p[frame]
    sec = time[frame]

    # Determine axis limit
    if isinstance(axis_lim_list, list):
        axis_lim = axis_lim_list[frame]  # Use frame-specific axis limit
    else:
        axis_lim = axis_lim_list  # Use single float value for a single frame

    # position and velocity vector of Sun
    #r_sun = p_t[2, 1:4]
    v_sun = p_t[2, 4:7]
    
    # position vector of Earth (as a proxy of HST)
    r_earth = p_t[3, 1:4]

    # calculate the two basis vectors of the projection plane
    l1 = np.cross(r_earth, v_sun)
    l2 = np.cross(l1, r_earth)
    l1 = l1 / np.linalg.norm(l1)
    l2 = l2 / np.linalg.norm(l2)

    # create an array to record coordinates of particles projected onto the plane
    p_projected = np.zeros((len(p_t),3), dtype=float)
    p_projected[:, 0] = p_t[:, 0]   # copy the column of particle ID   
    for i in range(len(p_t)):
        p_projected[i,1] = np.dot(l2, p_t[i, 1:4])
        p_projected[i,2] = np.dot(l1, p_t[i, 1:4])

    # Create a new figure for this frame
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # plot Didymos and Dimorphos
    ax.scatter(p_projected[0, 1], p_projected[0, 2], c='red', s=10, zorder=3, label='Didymos')
    ax.scatter(p_projected[1, 1], p_projected[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

    # plot dust particles
    ax.scatter(p_projected[4:, 1], p_projected[4:, 2], c='k', s=2)

    # Plot Sun direction relative to Didymos System Barycenter
    sun_x, sun_y = p_projected[2, 1], p_projected[2, 2]
    sun_distance = (sun_x**2 + sun_y**2) ** 0.5
    arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')

    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_xlabel('x / km')
    ax.set_ylabel('y / km')
    ax.grid()
    ax.legend()
   
    # Set title
    ax.set_title(f't = {sec/86400:.2f} days   Dust radius r = 1 mm')

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return frame_filename

def render_phase_angle_histogram(frame):
    """
    Renders a histogram of phase angles for a given frame and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.

    Returns:
        str: The path of the saved histogram image.
    """
    global shared_data
    data_p = shared_data['data_p']
    time = shared_data['time']
    output_dir = shared_data['output_dir']

    print(f"Rendering phase angle histogram for frame {frame}...", end="\n", flush=True)
    p_t = data_p[frame]
    sec = time[frame]

    # Get coordinates
    r_sun = p_t[2, 1:4]  # [x, y, z] of the Sun
    r_earth = p_t[3, 1:4]  # [x, y, z] of the Earth
    r_dust = p_t[4:, 1:4]  # [x, y, z] of all dust particles

    # Compute vectors
    vec_sun = r_sun - r_dust  # Vectors from Sun to dust particles
    vec_earth = r_earth - r_dust  # Vectors from Earth to dust particles

    # Compute norms
    norm_sun = np.linalg.norm(vec_sun, axis=1)  # Magnitudes of Sun vectors
    norm_earth = np.linalg.norm(vec_earth, axis=1)  # Magnitudes of Earth vectors

    # Compute dot products
    dot_product = np.einsum('ij,ij->i', vec_sun, vec_earth)  # Dot product of Sun and Earth vectors

    # Compute phase angles
    cos_phase_angle = dot_product / (norm_sun * norm_earth)  # Cosine of phase angle
    phase_angle = np.arccos(np.clip(cos_phase_angle, -1.0, 1.0))  # Phase angle in radians
    phase_angle_deg = np.degrees(phase_angle)
    
    # Compute histogram of phase angles
    counts, bin_edges = np.histogram(phase_angle_deg, bins=200, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Create distribution of phase angles
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(bin_centers, counts, color='blue', lw=2, label='')
    #plt.bar(bin_centers, counts, width=bin_edges[1]-bin_edges[0])
    ax.set_yscale('log')
    ax.set_xlabel('Phase Angle (degrees)')
    ax.set_ylabel('Number of Particles')
    ax.set_title(f'Phase Angle Distribution (t = {sec/86400:.2f} days)')
    
    # Save the histogram as a PNG image
    histogram_filename = os.path.join(output_dir, f"phase_angle_hist_{frame:04d}.png")
    plt.savefig(histogram_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return histogram_filename
