import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import map_coordinates
from coordinates import position_dict, day_hstfile_mapping

# --- 1. Data Loading and Preprocessing Function ---

def load_and_preprocess_fits(filepath, day_key):
	"""
	Loads FITS data, cleans it, applies log scaling, and calculates 
	the coordinate system parameters (pixel size, extent).

	Args:
		filepath (str): Full path to the FITS file.
		day_key (str): Key used to retrieve coordinate data from position_dict.

	Returns:
			tuple: (log10_data, extent, pixel_km, nx, ny)
	"""
	try:
		with fits.open(filepath) as hdul:
			hst_data = hdul[0].data
			header = hdul[0].header
	except FileNotFoundError:
		print(f"Error: FITS file not found at {filepath}")
		return None, None, None, None, None

	# Clean data, set background to 1e-20, and logarithmic brightness scale
	hst_data[hst_data < 0] = 0
	hst_data[np.isnan(hst_data)] = 0
	log10_hst = np.log10(hst_data + 1e-20)

	# --- Calculate Coordinate System ---
	
	# Coordinates are in 'Sun Body Center' - ensure day_key is correct
	r_Didy_sys_bary = position_dict['Didy_sys_bary'][day_key]
	r_Hubble = position_dict['Hubble'][day_key]
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

	# Determine the extent of the image for correct plotting
	x_min, x_max = x_km.min(), x_km.max()
	y_min, y_max = y_km.min(), y_km.max()
	extent = [x_min, x_max, y_min, y_max]
	
	return log10_hst, extent, pixel_km, nx, ny



def load_and_preprocess_h5(file_path, data_name):
	# 'r' means open the file for reading
	try:
		with h5py.File(file_path, 'r') as f:
			# Access the dataset by its name
			dset = f[data_name]
			# Read the entire dataset into a NumPy array in memory
			loaded_data = dset[:]
	except FileNotFoundError:
		print(f"Error: h5 file not found at {filepath}")
		return None, None, None

	ny, nx = loaded_data.shape
	log10_data = np.log10(loaded_data + 1e-20)

	return log10_data, nx, ny


# --- 2. Interactive Line Drawing and Sampling Function ---

# The sampling function remains largely the same, but it now expects the 
# data dimensions (Nx, Ny) and extent parameters separately if they were different. 
# Since we assume the system is identical, we pass the necessary parameters derived 
# from the first image.

def measure_tail_intensity(image_array_1, image_array_2, extent, nx, ny, output_dir):
	"""
	Displays image_array_1, prompts user to select two points, 
	samples intensity from both arrays, and plots results to compare.
	"""
	
	fig, ax = plt.subplots(figsize=(12, 12))

	# Display Image 1 for selection
	img_plot = ax.imshow(image_array_1, origin='lower', cmap='cividis',
	                    vmin=-8, vmax=np.nanmax(image_array_1), extent=extent)
	
	#cbar = fig.colorbar(img_plot, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
	#cbar.set_label('$\log_{10}$ (Pixel Value)')

	ax.set_xlabel('X Distance (km)')
	ax.set_ylabel('Y Distance (km)')
	ax.set_title('HST Image: Click to define tail line (Start and End)')

	print("--- INSTRUCTIONS ---")
	print("Please click two points on the image (Image 1): ")
	print("1. Start point of the tail axis.")
	print("2. End point of the tail axis.")
	print("--------------------")

	try:
		points = plt.ginput(2, timeout=90)
		plt.close(fig)
	except RuntimeError:
		print("\nSelection timed out or window closed.")
		plt.close(fig)
		return

	if len(points) < 2:
		print("\nError: Did not select two points.")
		return

	x1_km, y1_km = points[0]
	x2_km, y2_km = points[1]
	
	# --- Convert km coordinates to array indices (pixels) ---
	x_min, x_max, y_min, y_max = extent
	
	# Array indices (col, row)
	x1_pix = (x1_km - x_min) * (nx / (x_max - x_min))
	y1_pix = (y1_km - y_min) * (ny / (y_max - y_min))
	
	x2_pix = (x2_km - x_min) * (nx / (x_max - x_min))
	y2_pix = (y2_km - y_min) * (ny / (y_max - y_min))

	# --- Sample Pixel Values ---
	num_samples = 500
	rows = np.linspace(y1_pix, y2_pix, num_samples) 
	cols = np.linspace(x1_pix, x2_pix, num_samples)
	coords = [rows, cols]
	
	# Calculate the distance along the tail for the x-axis of the final plot
	distance_km = np.sqrt((np.linspace(x1_km, x2_km, num_samples) - x1_km)**2 + 
												(np.linspace(y1_km, y2_km, num_samples) - y1_km)**2)

	# Sample intensity profile from Image 1 and Image 2
	intensity_profile_1 = map_coordinates(image_array_1, coords, order=1)
	intensity_profile_2 = map_coordinates(image_array_2, coords, order=1)

	# --- 3. Plot the Intensity Profiles ---
	
	plt.figure(figsize=(10, 2))
	ftsize = 12

	plt.plot(distance_km, intensity_profile_1, 'k-', linewidth=2, label='HST observation')
	plt.plot(distance_km, intensity_profile_2, 'r--', linewidth=2, label='Model fit')
	#plt.plot(distance_km, intensity_profile_1, 'ko', markersize=2, label='HST observation')

	plt.xlabel('Distance Along Tail (km)', fontsize=ftsize+2)
	plt.ylabel('($\log_{10}$ Pixel Value)', fontsize=ftsize)
	plt.title(f'Intensity Profile Comparison Along Selected Tail Axis', fontsize=ftsize+4)
	plt.grid(True, linestyle='--', alpha=0.7)
	plt.legend(fontsize=ftsize)
	plt.savefig(os.path.join(output_dir, f'tail_intensity_profile.png'), dpi=300, bbox_inches='tight')
	plt.show()


# --- 4. Execution Block ---

if __name__ == "__main__":
	# Define file paths and coordinate key
	day_code = 'day_11.86' # Key for coordinate data
	FILEPATH_1 = os.path.join("/home/linfel/linfel_data/hst_raw_JianyangLi/", day_hstfile_mapping[day_code])
	FILEPATH_2 = f"/home/linfel/linfel_data/shortterm_anal/{day_code}_260-380/I_fit.h5"
	output_dir = f"/home/linfel/linfel_data/shortterm_anal/{day_code}_260-380"

	# Load and process Image 1 (This sets the primary coordinate system)
	log10_hst, extent_1, pixel_km_1, nx_1, ny_1 = load_and_preprocess_fits(FILEPATH_1, day_code)

	# Load and process Image 2
	log10_I_fit, nx_2, ny_2 = load_and_preprocess_h5(FILEPATH_2, 'intensity')

	if log10_hst is not None and log10_I_fit is not None:
		# Ensure both images are the same size before proceeding
		if log10_hst.shape != log10_I_fit.shape:
			print("Error: The two image arrays must have the same shape for comparison.")
		else:
			# Execute the main measurement and plotting function
			measure_tail_intensity(log10_hst, log10_I_fit, extent_1, nx_1, ny_1, output_dir)
