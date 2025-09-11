import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # Necessary for '3d' projection


def plot_particles_3d(filename, x_lim=None, y_lim=None, z_lim=None, point_size=5):
    """
    Reads particle data and creates an interactive 3D scatter plot.
    """
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Extract columns for position and velocity
    x_pos = data[:, 1]
    y_pos = data[:, 2]
    z_pos = data[:, 3]
    vx = data[:, 4]
    vy = data[:, 5]
    vz = data[:, 6]

    # Calculate the total speed for each particle
    #speed = np.sqrt(vx**2 + vy**2 + vz**2)

    # --- 3D Plotting Setup ---
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d') # <-- Create a 3D subplot

    # Create the 3D scatter plot
    scatter = ax.scatter(x_pos, y_pos, z_pos, c=vy, cmap='viridis', s=point_size,
		                     alpha=0.8, vmin=-10., vmax=2)

    # Add a color bar
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label('Speed (cm/s)', fontsize=12)

    # Set labels and title
    ax.set_xlabel('x (cm)', fontsize=12)
    ax.set_ylabel('y (cm)', fontsize=12)
    ax.set_zlabel('z (cm)', fontsize=12) # <-- Add Z label
    ax.set_title('3D Particle Positions Colored by Speed', fontsize=14)

    # Apply axis limits if they are provided
    if x_lim: ax.set_xlim(x_lim)
    if y_lim: ax.set_ylim(y_lim)
    if z_lim: ax.set_zlim(z_lim)

    plt.show()

def plot_vx(filename, x_lim=None, y_lim=None, point_size=10):
    """
    Reads particle data from a file and creates a scatter plot of x vs. y,
    with points colored by their total speed.

    Args:
        filename (str): Path to the data file.
        x_lim (tuple, optional): A tuple (xmin, xmax) for the x-axis limit.
        y_lim (tuple, optional): A tuple (ymin, ymax) for the y-axis limit.
        point_size (int, optional): The size of the markers in the scatter plot.
    """
    # Load the data from the text file, skipping the header if it exists
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Error reading file: {e}")
        print("Please ensure your file is a plain text file with numeric, space-separated values.")
        return

    # Extract columns for position and velocity
    # Fields: 0:ID, 1:x, 2:y, 3:z, 4:vx, 5:vy, 6:vz, 7:mass, 8:density
    x_pos = data[:, 1]
    y_pos = data[:, 2]
    vx = data[:, 4]
    vy = data[:, 5]
    vz = data[:, 6]

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    scatter = ax.scatter(x_pos, y_pos, c=vx, cmap='viridis', s=point_size,
		                     alpha=0.8, vmin=-10., vmax=10.)

    # Add a color bar to show the speed scale
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r'$v_x$ (cm/s)', fontsize=12)

    # Set labels and title
    ax.set_xlabel('x (cm)', fontsize=12)
    ax.set_ylabel('y (cm)', fontsize=12)
    ax.set_title(r'Particle Positions Colored by $v_x$', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_aspect('equal', adjustable='box')

    # Apply axis limits if they are provided
    if x_lim:
        ax.set_xlim(x_lim)
    if y_lim:
        ax.set_ylim(y_lim)

    plt.show()


def plot_vy(filename, x_lim=None, y_lim=None, point_size=10):
    """
    Reads particle data from a file and creates a scatter plot of x vs. y,
    with points colored by their total speed.

    Args:
        filename (str): Path to the data file.
        x_lim (tuple, optional): A tuple (xmin, xmax) for the x-axis limit.
        y_lim (tuple, optional): A tuple (ymin, ymax) for the y-axis limit.
        point_size (int, optional): The size of the markers in the scatter plot.
    """
    # Load the data from the text file, skipping the header if it exists
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        print(f"Error reading file: {e}")
        print("Please ensure your file is a plain text file with numeric, space-separated values.")
        return

    # Extract columns for position and velocity
    # Fields: 0:ID, 1:x, 2:y, 3:z, 4:vx, 5:vy, 6:vz, 7:mass, 8:density
    x_pos = data[:, 1]
    y_pos = data[:, 2]
    vx = data[:, 4]
    vy = data[:, 5]
    vz = data[:, 6]

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    scatter = ax.scatter(x_pos, y_pos, c=vy, cmap='viridis', s=point_size,
		                     alpha=0.8, vmin=-10., vmax=0.)

    # Add a color bar to show the speed scale
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r'$v_y$ (cm/s)', fontsize=12)

    # Set labels and title
    ax.set_xlabel('x (cm)', fontsize=12)
    ax.set_ylabel('y (cm)', fontsize=12)
    ax.set_title(r'Particle Positions Colored by $v_y$', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_aspect('equal', adjustable='box')

    # Apply axis limits if they are provided
    if x_lim:
        ax.set_xlim(x_lim)
    if y_lim:
        ax.set_ylim(y_lim)

    plt.show()


if __name__ == '__main__':

    # --- USAGE EXAMPLES ---
    print("\nDisplaying zoomed-in plot...")
#    plot_particles_3d("/nuke/linfel/Ejecta/data_0Pa_160s_dimor_removed.txt", 
#                      x_lim=(-200e2, 200e2), y_lim=(-200e2, 200e2), z_lim=(-200e2, 200e2))
    
    plot_vx("/nuke/linfel/Ejecta/data_0Pa_160s_dimor_removed.txt", 
                      x_lim=(-200e2, 200e2), y_lim=(-200e2, 200e2))
