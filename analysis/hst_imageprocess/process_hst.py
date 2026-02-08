import os
import numpy as np
import matplotlib.pyplot as plt
import json
from matplotlib.path import Path

def select_area_and_get_indices(image_array, extent, nx, ny, mask_filepath=None):
    """
    1. Prompts the user to draw a polygon on the displayed image (if mask not loaded).
    2. Calculates and returns a 1D boolean mask for pixels inside the polygon.
    3. Optionally saves/loads the boolean mask to/from a .npy file.

    Args:
        image_array (np.ndarray): The 2D image data array (e.g., log10_hst).
        extent (list): [xmin, xmax, ymin, ymax] in km.
        nx (int): Number of columns (x-dimension size) in the array.
        ny (int): Number of rows (y-dimension size) in the array.
        mask_filepath (str, optional): Path to save/load the 1D boolean mask.

    Returns:
        np.ndarray: A 1D boolean mask (size nx*ny) where True indicates the pixel is inside the polygon.
        list: The list of (x_km, y_km) vertices used to create the mask (for plotting).
    """
    
    # --- Load Mask from File ---
    if mask_filepath and os.path.exists(mask_filepath):
        try:
            # Load the 1D boolean mask
            polygon_mask_1d = np.load(mask_filepath)
            print(f"Loaded 1D boolean mask from {mask_filepath}")
            # We don't have the vertices, so we return an empty list for plotting
            return polygon_mask_1d, [] 
        except Exception as e:
            print(f"Error loading mask file: {e}. Proceeding with interactive selection.")

    # --- Interactive Drawing (If mask not loaded) ---
    # (The plotting and ginput logic remains the same as before to get vertices_km)
    
    # ... [Interactive drawing logic to get vertices_km, requires figure setup] ...
    # Placeholder for brevity, assume vertices_km is obtained interactively:
    
    # --- Start Interactive Drawing ---
    fig, ax = plt.subplots(figsize=(10, 10))
    # Display the original (not log10) data for visual selection if possible, or use log10_hst
    ax.imshow(image_array, origin='lower', cmap='cividis', vmin=-14, vmax=np.nanmax(image_array), extent=extent)
    ax.set_title('Click to define polygon vertices. Press Enter when finished.')
    plt.show(block=False) 
    
    print("Click to define vertices. Press Enter (or double-click) to complete the polygon.")
    try:
        vertices_km = plt.ginput(n=-1, timeout=180, mouse_add=1, mouse_pop=3, mouse_stop=2)
        plt.close(fig)
    except RuntimeError:
        print("Selection timed out or window closed.")
        plt.close(fig)
        return np.array([]), []

    if len(vertices_km) < 3:
        print("Requires at least 3 points to form a polygon.")
        return np.array([]), []
    # --- End Interactive Drawing ---

    x_min, x_max, y_min, y_max = extent
    
    # --- 2. Convert km vertices to pixel indices ---
    x_range, y_range = x_max - x_min, y_max - y_min
    vertices_pix = []
    for x_km, y_km in vertices_km:
        x_pix = (x_km - x_min) * (nx / x_range)
        y_pix = (y_km - y_min) * (ny / y_range)
        vertices_pix.append((x_pix, y_pix))

    # --- 3. Identify Pixels Inside the Polygon and Flatten ---
    y_indices, x_indices = np.mgrid[0:ny, 0:nx]
    all_pixels = np.vstack((x_indices.ravel(), y_indices.ravel())).T # (N*M, 2) array of (col, row)
    
    polygon_path = Path(vertices_pix)
    polygon_mask_1d = polygon_path.contains_points(all_pixels) # This is the 1D boolean mask

    # --- 4. Save the Mask ---
    if mask_filepath:
        np.save(mask_filepath, polygon_mask_1d)
        print(f"1D boolean mask saved to {mask_filepath}")

    return polygon_mask_1d, vertices_km


def plot_hst_image_with_selection(log10_hst, extent, selected_vertices,  output_dir, day_code):
    """
    Plots the HST log-scaled image, overlays the selected polygon area,
    and saves the resulting figure.
    """
    ny, nx = log10_hst.shape
    x_min, x_max, y_min, y_max = extent

    # Create coordinate arrays matching your original logic
    x_km_edges = np.linspace(x_min, x_max, nx + 1)
    y_km_edges = np.linspace(y_min, y_max, ny + 1)
    X, Y = np.meshgrid(x_km_edges, y_km_edges)

    # 2. Plot HST observation
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Use pcolormesh with the calculated grid
    pc = ax.pcolormesh(X, Y, log10_hst, cmap='cividis', shading='auto', vmin=-10, vmax=np.nanmax(log10_hst))

    # 3. Plot the Polygon Vertices
    if selected_vertices:
        # Separate X and Y coordinates of the vertices
        x_verts = [v[0] for v in selected_vertices]
        y_verts = [v[1] for v in selected_vertices]

        # Close the loop for plotting the polygon
        x_verts.append(x_verts[0])
        y_verts.append(y_verts[0])

        # Overlay the polygon outline
        ax.plot(
            x_verts,
            y_verts,
            color='red',
            linewidth=2,
            linestyle='-',
            label='Selected Area',
            zorder=10 # Ensure the polygon is visible on top
        )
        
        # Plot the individual vertices as markers
        ax.plot(
            x_verts[:-1],
            y_verts[:-1],
            'ro', # Red circles
            markersize=5,
            zorder=11
        )

    # 4. Final Plot Formatting
    cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.8, aspect=30)
    cbar.set_label(r'$\log_{10}$(Pixel Value)')

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('X Distance (km)')
    ax.set_ylabel('Y Distance (km)')
    ax.set_title(f'HST Image (Selected Area Highlighted)')
    
    # You might want to remove the grid if the polygon lines are enough visual guidance
    # ax.grid(True, linestyle='--', linewidth=0.1, color='white', alpha=0.7)
    
    plt.legend()
    plt.tight_layout()

    # 5. Save and Close
    filename = os.path.join(output_dir, f'hst_image_{day_code}_with_selection.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {filename}")
    plt.close()


def plot_selected_region(log10_hst, polygon_mask_1d, extent, output_dir, day_code):
    """
    Plots only the selected region of the HST image using the boolean mask,
    against a transparent background.

    Args:
        log10_hst (np.ndarray): 2D array of log10(Pixel Value) (shape: ny, nx).
        polygon_mask_1d (np.ndarray): 1D boolean mask (size ny*nx) for the selected area.
        extent (list): [xmin, xmax, ymin, ymax] in km for image scaling.
        utc_mid (str): Date/time string for the plot title.
        output_dir (str): Directory to save the output file.
        day_code (str): Code for the filename.
    """
    ny, nx = log10_hst.shape
    
    # 1. Reshape the 1D mask back to 2D
    # The mask is True for selected pixels.
    polygon_mask_2d = polygon_mask_1d.reshape(ny, nx)

    # 2. Create the masked data array
    # We create a copy of the original data.
    masked_data = log10_hst.copy()
    
    # The inverse mask (logical NOT) selects the *unwanted* pixels.
    inverse_mask = ~polygon_mask_2d
    
    # Set all unwanted pixels to NaN (which makes them transparent in Matplotlib)
    masked_data[inverse_mask] = np.nan
    
    # 3. Setup Plot
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Create coordinate arrays for pcolormesh edges
    x_min, x_max, y_min, y_max = extent
    x_km_edges = np.linspace(x_min, x_max, nx + 1)
    y_km_edges = np.linspace(y_min, y_max, ny + 1)
    X, Y = np.meshgrid(x_km_edges, y_km_edges)

    # 4. Plot the Masked Data
    # Use the original vmin/vmax for consistent color scale
    pc = ax.pcolormesh(
        X, Y, masked_data, 
        cmap='cividis', 
        shading='auto', 
        vmin=-14, # Using the typical dim limit from your previous code
        vmax=np.nanmax(log10_hst) # Max value from the original data
    )

    # 5. Final Plot Formatting
    cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.8, aspect=30)
    cbar.set_label(r'$\log_{10}$(Pixel Value)')

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('X Distance (km)')
    ax.set_ylabel('Y Distance (km)')
    ax.set_title(f'HST Image Selection on ')
    
    # Set the background color to be explicit (e.g., white or black)
    ax.set_facecolor('black') 
    
    plt.tight_layout()

    # 6. Save and Close
    filename = os.path.join(output_dir, f'hst_image_{day_code}_selected_region.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Plot of selected region saved to: {filename}")
    plt.close()
