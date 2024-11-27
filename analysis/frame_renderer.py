import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os

# Define global variables for shared data
shared_data = {}

def init_worker(data_p, time, axis_lim_list, output_dir):
    """
    Initialize global shared data for each worker process.
    This function is called once per worker during pool initialization.
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
    Renders a single frame and saves it as a PNG file.

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
    ax.set_title(f't = {sec/86400:.2f} days   Dust particle radius r = 1 mm')

    # Save the frame as a PNG image
    frame_filename = os.path.join(output_dir, f"frame_{frame:04d}.png")
    plt.savefig(frame_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return frame_filename
