import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from coordinates import position_dict
from scipy.ndimage import map_coordinates

# --- 1. Load the Image Data (Placeholder) ---
# Replace this section with the code that loads your 2000x2000 pixel data array.
# The image you provided is a log-scaled plot, so the array should contain the 
# raw log10(Pixel Value) data, or the original intensity values.

# For this example, we'll create a fake 2000x2000 array to simulate the image data.
# ASSUMPTION: 'image_data' is a 2D NumPy array (2000, 2000) containing 
# the log10(Pixel Value) data displayed in your plot.

# --- Replace this placeholder with your actual image data loading ---
# Example: image_data = np.load('your_hst_data.npy')
# OR if you are loading the image file itself and need to extract the array:
# from PIL import Image
# img = Image.open('hst_image_day_64.44.jpg').convert('L') 
# image_data = np.array(img.getdata()).reshape(img.size) # This only gets visual data.
# You MUST use the astronomical data array (FITS/HDF5) that generated the plot.
# Load FITS image
with fits.open('/home/linfel/linfel_data/hst_raw_JianyangLi/16674/stack_18_long.fits.fits') as hdul:
    hst_data = hdul[0].data  # assuming image is in HDU 0
    header = hdul[0].header

# Extract metadata
orientat = header.get('ORIENTAT', 'N/A')
utc_mid = header.get('UTC-MID', 'N/A')

# clean data, set background to 1e-20, and logarithmic brightness scale
hst_data[hst_data < 0] = 0
hst_data[np.isnan(hst_data)] = 0
log10_hst = np.log10(hst_data + 1e-20)

# coordinates are in 'Sun Body Center'
r_Didy_sys_bary = position_dict['Didy_sys_bary']['day_5.70']
r_Hubble = position_dict['Hubble']['day_5.70']
target_distance = np.linalg.norm(r_Didy_sys_bary - r_Hubble)  # km

# Hubble pixel size
pixel_arcsec = 0.04
range_km = target_distance
pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad
pixel_km = 2 * range_km * np.tan(pixel_fov/2)

# Axes in km
ny, nx = hst_data.shape
x_km = np.arange(nx+1) * pixel_km
y_km = np.arange(ny+1) * pixel_km
x_km = x_km - x_km[-1]/2
y_km = y_km - y_km[-1]/2
# -------------------------------------------------------------------

# Placeholder: Create a noisy array with a simulated dust tail (for testing)
#IMAGE_SIZE = 100
#image_data = np.random.normal(loc=-13, scale=0.5, size=(IMAGE_SIZE, IMAGE_SIZE))
## Add a bright 'tail' diagonally
#y, x = np.mgrid[0:IMAGE_SIZE, 0:IMAGE_SIZE]
#center = IMAGE_SIZE // 2
#tail_brightness = np.exp(-((x - center) * 0.005 + (y - center) * 0.005)**2 / 0.1) * 10 
#image_data += tail_brightness 
# -------------------------------------------------------------------

# Determine the extent of the image for correct plotting (from your axes)
# Your image covers roughly +/- 5000 km, so the pixel size is 10000 km / 2000 pixels = 5 km/pixel
x_min, x_max = x_km.min(), x_km.max()
y_min, y_max = y_km.min(), y_km.max()
extent = [x_min, x_max, y_min, y_max]

# --- 2. Interactive Line Drawing and Sampling ---

def measure_tail_intensity(image_array, extent):
    """
    Displays the image, prompts user to select two points, samples intensity, and plots results.
    
    Args:
        image_array (np.ndarray): The 2D image data array (e.g., log10(Pixel Value)).
        extent (list): [xmin, xmax, ymin, ymax] for correct axes scaling.
    """
    
    fig, ax = plt.subplots(figsize=(10, 10))

    # Display the image using the defined extent
    img_plot = ax.imshow(image_array, origin='lower', cmap='gray', vmin=-14, vmax=-4, extent=extent)
    
    # Add colorbar
    cbar = fig.colorbar(img_plot, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
    cbar.set_label('$\log_{10}$ (Pixel Value)')

    # Set axes labels based on the plot context
    ax.set_xlabel('X Distance (km)')
    ax.set_ylabel('Y Distance (km)')
    ax.set_title('HST Image: Click to define tail line (Start and End)')

    print("--- INSTRUCTIONS ---")
    print("Please click two points on the image: ")
    print("1. Start point of the tail axis.")
    print("2. End point of the tail axis.")
    print("--------------------")

    # Get user input for two points (returns list of tuples: [(x1, y1), (x2, y2)])
    # The coordinates returned are in the plot's data units (km)
    try:
        points = plt.ginput(2, timeout=30)
        plt.close(fig) # Close the image window after selection
    except RuntimeError:
        print("\nSelection timed out or window closed.")
        plt.close(fig)
        return

    if len(points) < 2:
        print("\nError: Did not select two points.")
        return

    # Extract coordinates in kilometers
    x1_km, y1_km = points[0]
    x2_km, y2_km = points[1]
    print(f"\nLine defined from: ({x1_km:.1f}, {y1_km:.1f}) to ({x2_km:.1f}, {y2_km:.1f}) km")

    # --- 3. Convert km coordinates to array indices (pixels) ---
    
    # Calculate conversion factors
    # Assuming the array index 0 maps to x_min and index N-1 maps to x_max
    Nx, Ny = image_array.shape
    
    # Array indices (col, row)
    # The formula is: index = (km_value - min_km) * (N_pixels / total_km_range)
    x1_pix = (x1_km - x_min) * (Nx / (x_max - x_min))
    y1_pix = (y1_km - y_min) * (Ny / (y_max - y_min))
    
    x2_pix = (x2_km - x_min) * (Nx / (x_max - x_min))
    y2_pix = (y2_km - y_min) * (Ny / (y_max - y_min))

    # --- 4. Sample Pixel Values ---
    
    # Create the coordinates along the line
    num_samples = 500  # Number of points to sample along the line
    
    # Linearly space the row and column coordinates (indices)
    rows = np.linspace(y1_pix, y2_pix, num_samples) 
    cols = np.linspace(x1_pix, x2_pix, num_samples)
    
    # Calculate the distance along the tail for the x-axis of the final plot
    # Distance in km: Euclidean distance from the start point
    distance_km = np.sqrt((np.linspace(x1_km, x2_km, num_samples) - x1_km)**2 + 
                          (np.linspace(y1_km, y2_km, num_samples) - y1_km)**2)

    # Use map_coordinates to efficiently sample the image array
    # The coordinates for map_coordinates must be in (row, col) format
    # The order parameter '1' means linear interpolation (better than nearest neighbor)
    intensity_profile = map_coordinates(image_array, [rows, cols], order=1)

    # --- 5. Plot the Intensity Profile ---
    
    plt.figure(figsize=(10, 6))
    
    # Plot the sampled log10(Pixel Value) against distance along the tail
    plt.plot(np.log10(distance_km), intensity_profile, 'r-', linewidth=2)
    plt.xlabel('Log10 Distance Along Tail (km)')
    plt.ylabel('Measured Intensity ($\log_{10}$ Pixel Value)')
    plt.title(f'Intensity Profile Along Selected Tail Axis')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.show()

# Execute the function
measure_tail_intensity(log10_hst, extent)
