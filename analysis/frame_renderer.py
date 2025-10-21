import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LogNorm
import os
import numpy as np

# Define global variables for shared data
shared_data = {}

def init_worker(data_p, time, axis_lim_list, output_dir, group_labels=None):
    """
    Initialize global shared data for each worker process.
    This function is called once per worker during pool initialization.
    
    Args:
        data_p (list): A list of ndarray containing particle data for all frames.
        time (ndarray): Time data for all frames.
        axis_lim_list (list or float): List of axis limits for each frame (if rendering video) 
                                       or a single float value for one frame
                                       or None for default axis limit.
        output_dir (str): Directory to save the output image.
    """
    global shared_data
    shared_data['data_p'] = data_p
    shared_data['time'] = time
    shared_data['axis_lim_list'] = axis_lim_list
    shared_data['output_dir'] = output_dir
    shared_data['group_labels'] = group_labels

# Convert length unit for axis labels
def m_to_km(x, _):
    return f'{x / 1e3:.1f}'

def m_to_au(x, _):
    return f'{x / 1.495978707e11:.1f}'

def render_frame_topview(frame, dpi=100):
    """
    Renders a single frame of top view and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.

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
    
    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_aspect('equal', adjustable='box')
    
    # Plot Didymos and Dimorphos
    ax.scatter(p_t[0, 1], p_t[0, 2], c='red', s=10, zorder=3, label='Didymos')
    ax.scatter(p_t[1, 1], p_t[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

    # Plot dust particles
    speed = (p_t[4:, 4]**2 + p_t[4:, 5]**2 + p_t[4:, 6]**2)**0.5
    dust_sc = ax.scatter(p_t[4:, 1], p_t[4:, 2], c=speed, s=2, cmap='viridis', norm=LogNorm())
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

    ax.set_xlabel('x / km')
    ax.set_ylabel('y / km')
    ax.grid()
    ax.legend(loc='upper right')

    # Set title
    ax.set_title(f't = {sec/86400:.2f} days   Top View')

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"topview_frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return frame_filename

def render_frame_sideview(frame, dpi=100):
    """
    Renders a single frame of side view and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.

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

    # Create a new figure for this frame
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_aspect('equal', adjustable='box')
    
    # plot Didymos and Dimorphos
    ax.scatter(p_projected[0, 1], p_projected[0, 2], c='red', s=10, zorder=3, label='Didymos')
    ax.scatter(p_projected[1, 1], p_projected[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

    # Plot dust particles
    speed = (p_t[4:, 4]**2 + p_t[4:, 5]**2 + p_t[4:, 6]**2)**0.5
    dust_sc = ax.scatter(p_projected[4:, 1], p_projected[4:, 2], c=speed, s=2, cmap='viridis', norm=LogNorm())
    cb = fig.colorbar(dust_sc, ax=ax, shrink=0.8, aspect=20)
    cb.set_label(f'speed (m/s)', fontsize=12)
    
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
    ax.set_title(f't = {sec/86400:.2f} days   Side view positions')
    
    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"sideview_frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return frame_filename


def render_frame_HSTview(frame, dpi=100):
    """
    Renders a single frame of HST view and saves it as a PNG file.

    Args:
        frame (int): The frame index to render.

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

    # matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
    SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                              [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                              [-0.104839674791979,  0.124915784491013,  0.986612735258626]])
    # calculate sky north vector in 'Sun Body Center' and 
    # convert it to 'Didymos System Barycenter' frame
    obliq_earth = np.deg2rad(23.4392911)
    sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
    r_sky_north = SBC_rotate_DSB @ sky_north
    
    # position vector of Earth (as a proxy of HST)
    r_earth = p_t[3, 1:4]

    # calculate the two basis vectors of the projection plane
    l1 = np.cross(r_sky_north, r_earth)
    l2 = np.cross(r_earth, l1)
    l1 = l1 / np.linalg.norm(l1)
    l2 = l2 / np.linalg.norm(l2)

    # create an array to record coordinates of particles projected onto the plane
    p_projected = np.zeros((len(p_t),3), dtype=float)
    p_projected[:, 0] = p_t[:, 0]   # copy the column of particle ID   
    p_projected[:, 1] = np.dot(p_t[:, 1:4], l1)    # x coordinate
    p_projected[:, 2] = np.dot(p_t[:, 1:4], l2)    # y coordinate

    # Create a new figure for this frame
    fig, ax = plt.subplots(figsize=(8, 8))

    # Determine axis limit
    if axis_lim_list is None:
        print("No axis limits specified, using default limits.")
    elif isinstance(axis_lim_list, list):
        axis_lim = axis_lim_list[frame]  # Use frame-specific axis limit
        ax.set_xlim(-axis_lim, axis_lim)
        ax.set_ylim(-axis_lim, axis_lim)
    else:
        axis_lim = axis_lim_list         # Use single float value for a single frame  
        ax.set_xlim(-axis_lim, axis_lim)
        ax.set_ylim(-axis_lim, axis_lim)
    ax.set_aspect('equal', adjustable='box')
    
    # plot Didymos and Dimorphos
    ax.scatter(p_projected[0, 1], p_projected[0, 2], c='red', s=10, zorder=3, label='Didymos')
    ax.scatter(p_projected[1, 1], p_projected[1, 2], c='blue', s=8, zorder=3, label='Dimorphos')

    # plot dust particles
    ax.scatter(p_projected[4:, 1], p_projected[4:, 2], c='k', s=0.5, alpha=0.05)

    # Plot Sun direction relative to Didymos System Barycenter
    sun_x, sun_y = p_projected[2, 1], p_projected[2, 2]
    sun_distance = (sun_x**2 + sun_y**2) ** 0.5
    arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')

    # Heliocentric velocity direction of Didymos
    v_didy_heliocentric = p_t[0, 4:7] - p_t[2, 4:7]
    hvdd_x_proj = np.dot(v_didy_heliocentric, l1)
    hvdd_y_proj = np.dot(v_didy_heliocentric, l2)
    ratio = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1 / (hvdd_x_proj**2 + hvdd_y_proj**2) ** 0.5
    arrow_x = hvdd_x_proj * ratio
    arrow_y = hvdd_y_proj * ratio
    ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='cyan',
              label='Heliocentric velocity direction of Didymos')
    
    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlabel('x / km')
    ax.set_ylabel('y / km')
    ax.grid()
    ax.legend(loc='upper right')
   
    # Set title
    ax.set_title(f't = {sec/86400:.2f} days   HST Perspective')

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"hstview_frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return frame_filename

def render_single_HSTview_colorgroups(frame, alpha):
    """
    Renders a single frame of HST view, colorcoding particles of different groups.

    Args:
        frame (int): The frame index to render.

    Returns:
        str: The path of the saved frame image.
    """
    global shared_data
    data_p = shared_data['data_p']
    time = shared_data['time']
    axis_lim_list = shared_data['axis_lim_list']
    output_dir = shared_data['output_dir']
    group_labels = shared_data['group_labels']

    print(f"Rendering frame {frame}...", end="\r", flush=True)
    p_t = data_p[frame]
    sec = time[frame]

    # matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
    SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                              [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                              [-0.104839674791979,  0.124915784491013,  0.986612735258626]])
    # calculate sky north vector in 'Sun Body Center' and 
    # convert it to 'Didymos System Barycenter' frame
    obliq_earth = np.deg2rad(23.4392911)
    sky_north = np.array([0, np.sin(obliq_earth), np.cos(obliq_earth)])
    r_sky_north = SBC_rotate_DSB @ sky_north

    # position vector of Earth (as a proxy of HST)
    r_earth = p_t[3, 1:4]

    # calculate the two basis vectors of the projection plane
    l1 = np.cross(r_sky_north, r_earth)
    l2 = np.cross(r_earth, l1)
    l1 = l1 / np.linalg.norm(l1)
    l2 = l2 / np.linalg.norm(l2)

    # create an array to record coordinates of particles projected onto the plane
    p_projected = np.zeros((len(p_t),4), dtype=float)
    p_projected[:, 0] = p_t[:, 0]   # copy the column of particle ID
    p_projected[:, -1] = p_t[:, -1] # copy the column of particle group      
    p_projected[:, 1] = np.dot(p_t[:, 1:4], l1)    # x coordinate
    p_projected[:, 2] = np.dot(p_t[:, 1:4], l2)    # y coordinate
    
    # Create a new figure for this frame
    fig, ax = plt.subplots(figsize=(8,8))
    font = 15

    # Determine axis limit
    if axis_lim_list is None:
        print("No axis limits specified, using default limits.")
    elif isinstance(axis_lim_list, list):
        if len(axis_lim_list)==4:         # set upper, bottom, left, right axis limits
            ax.set_xlim(axis_lim_list[0], axis_lim_list[1])
            ax.set_ylim(axis_lim_list[2], axis_lim_list[3])
            ax.set_aspect('equal', adjustable='box')
    elif isinstance(axis_lim_list, float):
        axis_lim = axis_lim_list         # Use single float value for a single frame  
        ax.set_xlim(-axis_lim, axis_lim)
        ax.set_ylim(-axis_lim, axis_lim)
        ax.set_aspect('equal', adjustable='box')
    
    # plot Didymos and Dimorphos
    sc_didy = ax.scatter(p_projected[0, 1], p_projected[0, 2], c='k', s=10, zorder=3, label='Didymos')
    sc_dimor = ax.scatter(p_projected[1, 1], p_projected[1, 2], c='blue', s=8, zorder=2, label='Dimorphos')

    # plot dust particles in groups
    for group, charac in group_labels.items():
        group_indices = p_projected[:, -1] == group  # group info is in the last column
        ax.scatter(p_projected[group_indices, 1],
                   p_projected[group_indices, 2],
                   c=charac[1],
                   s=0.5,
                   label=f'{charac[0]}',
                   alpha=alpha)

    # Plot Sun direction relative to Didymos System Barycenter
    sun_x, sun_y = p_projected[2, 1], p_projected[2, 2]
    sun_distance = (sun_x**2 + sun_y**2) ** 0.5
    arrow_x = sun_x / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    arrow_y = sun_y / sun_distance * (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.1
    qv_sun = ax.quiver(0, 0, arrow_x, arrow_y, angles='xy', scale_units='xy', scale=1, width=0.005, color='orange', label='Sun Direction')

    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_km))
    ax.set_xlabel('x / km',fontsize=font)
    ax.set_ylabel('y / km',fontsize=font)
    ax.grid()

    # configure legend
    handles = [           # create new handles for dusts
        plt.Line2D([0], [0], marker='o', color=charac[1], markersize=2, linestyle='None', label=f'{charac[0]}')
        for group, charac in group_labels.items()]
    handles.extend([sc_didy, sc_dimor, qv_sun])
    ax.legend(handles=handles, loc='best', fontsize=font-3, ncol=4)
   
    # Set title
    ax.set_title(f't = {sec/86400:.2f} days', size=font)

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"day_{sec/86400:.3f}.png")
    plt.savefig(frame_filename, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close(fig)
    return frame_filename

def render_location_sun_center(frame):
    """
    Renders a single frame of locations of Sun, Earth, Didymos system, velocity direction displayed.

    Args:
        frame (int): The frame index to render.

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
    assert isinstance(axis_lim_list, float), "Data type error. 'axis_lim_list' is not a float."
    axis_lim = axis_lim_list

    # position and velocity vector of Sun and Earth (as a proxy of HST)
    r_sun = p_t[2, 1:4]
    v_sun = p_t[2, 4:7]
    r_earth = p_t[3, 1:4]
    v_earth = p_t[3, 4:7]

    # matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
    SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                              [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                              [-0.104839674791979,  0.124915784491013,  0.986612735258626]])
    DSB_rotate_SBC = np.linalg.inv(SBC_rotate_DSB)
    
    # calculate position and velocity vectors of DSB and Earth in SBC reference frame
    r_DSB_1= -np.dot(DSB_rotate_SBC, r_sun.T)
    v_DSB_1= -np.dot(DSB_rotate_SBC, v_sun.T)
    r_earth_1 = np.dot(DSB_rotate_SBC, r_earth.T) + r_DSB_1
    v_earth_1 = np.dot(DSB_rotate_SBC, v_earth.T) + v_DSB_1
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # plot positions of Sun, Earth, DSB
    ax.scatter(0, 0, c='red', s=15, zorder=3, label='Sun')
    ax.scatter(r_earth_1[0], r_earth_1[1], c='green', s=10, zorder=3, label='Earth')
    ax.scatter(r_DSB_1[0], r_DSB_1[1], c='k', s=8, zorder=3, label='Didymos system')

    # Plot DSB velocity direction
    v_arrow = v_DSB_1 / np.linalg.norm(v_DSB_1) * axis_lim * 0.1
    ax.quiver(r_DSB_1[0], r_DSB_1[1], v_arrow[0], v_arrow[1], angles='xy', scale_units='xy',
              scale=1, width=0.002, color='grey')

    # Plot Earth velocity direction
    v_arrow = v_earth_1 / np.linalg.norm(v_earth_1) * axis_lim * 0.1
    ax.quiver(r_earth_1[0], r_earth_1[1], v_arrow[0], v_arrow[1], angles='xy', scale_units='xy',
              scale=1, width=0.002, color='grey')

    # Customize axis
    ax.xaxis.set_major_formatter(FuncFormatter(m_to_au))
    ax.yaxis.set_major_formatter(FuncFormatter(m_to_au))
    ax.set_xlim(-axis_lim, axis_lim)
    ax.set_ylim(-axis_lim, axis_lim)
    ax.set_xlabel('x / au')
    ax.set_ylabel('y / au')
    ax.grid()
    ax.legend(loc='upper right')
   
    # Set title
    ax.set_title(f't = {sec/86400:.2f} days')

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
