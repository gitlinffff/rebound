from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

"""
HST
"""
# Load FITS image
with fits.open('/home/linfel/linfel_data/hst_raw_JianyangLi/17289/stack_all_long.fits') as hdul:
    hst_data = hdul[0].data  # assuming image is in HDU 0
    header = hdul[0].header

# Extract metadata
orientat = header.get('ORIENTAT', 'N/A')
utc_mid = header.get('UTC-MID', 'N/A')

# the time of the hst image
hst_time = datetime.strptime(utc_mid, "%Y-%m-%dT%H:%M:%S.%f")

# Time of the impact
impact_datetime = "2022-09-26T23:14:24.1830"
impact_time = datetime.strptime(impact_datetime, "%Y-%m-%dT%H:%M:%S.%f")

# calculate the time difference
time_difference = hst_time - impact_time
days_after_impact = time_difference.total_seconds() / 86400.
print(f"Time of HST data: T0+{days_after_impact:.2f} days")

# clean data, set background to 1e-20, and logarithmic brightness scale
hst_data[hst_data < 0] = 0
hst_data[np.isnan(hst_data)] = 0
log10_hst = np.log10(hst_data + 1e-20)

# coordinates are in 'Sun Body Center'
r_Didy_sys_bary = np.array([-1.252938965664409E+08, 1.736400744097055E+08, 1.017805938869286E+07])  # km
r_Hubble = np.array([-1.058142769367994E+08, 1.027338649631926E+08, -2.109203366145492E+03])        # km
target_distance = np.linalg.norm(r_Didy_sys_bary - r_Hubble)  # km

# Hubble pixel size
pixel_arcsec = 0.04
range_km = target_distance
pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad
pixel_km = 2 * range_km * np.tan(pixel_fov/2)

# Axes limits in km
ny, nx = hst_data.shape
x_km = np.arange(nx+1) * pixel_km
y_km = np.arange(ny+1) * pixel_km
x_km = x_km - x_km[-1]/2
y_km = y_km - y_km[-1]/2

"""
Prepare Simulation Pixel Data
"""
import os, pickle
from ReadParticle import read_particle_frames
import miepython as mie

def rotate_coords(x, y, angle_deg, center=(0, 0)):
    """
    Rotate coordinates (x, y) by angle_deg around a given center.

    Parameters:
    - x, y: float or np.array of coordinates
    - angle_deg: rotation angle in degrees (positive = counterclockwise)
    - center: (x0, y0) tuple for rotation center (default is origin)

    Returns:
    - x_rot, y_rot: rotated coordinates
    """
    angle_rad = np.radians(angle_deg)

    x0, y0 = center

    x_rot = np.cos(angle_rad) * (x - x0) - np.sin(angle_rad) * (y - y0) + x0
    y_rot = np.sin(angle_rad) * (x - x0) + np.cos(angle_rad) * (y - y0) + y0

    return x_rot, y_rot

# data path
data_rootdir = "/home/linfel/linfel_data/ejecta_exp_datahigh/day83.77_pkl"
filenames = [os.path.join(data_rootdir, f"particle_{i:03d}_day83.77.pkl") for i in range(9, 22)]

# Ensure output directory exists for saving frames
output_dir = "/home/linfel/linfel_data/ejecta_exp_datahigh/postprocess/weight_fit_result"
os.makedirs(output_dir, exist_ok=True)

# matrix convert vector from 'Sun Body Center' to 'Didymos System Barycenter'
SBC_rotate_DSB = np.array([[-0.703595792353257, -0.710438191316344, 0.015183454875444],
                          [ 0.702824020343859, -0.692584740738747,  0.16237232930379 ],
                          [-0.104839674791979,  0.124915784491013,  0.986612735258626]])

# sky north vector and convert it to 'Didymos System Barycenter' frame
sky_north = np.array([0, 0, 1])
r_sky_north = np.dot(SBC_rotate_DSB, sky_north.T)

# optical parameters
m = 1.7 - 0.01j           # refractive index of particle
lambda0 = 500e-9  # wavelength in vacuum (m)

# x and y limits of the intensity map 
xedges = x_km * 1e3
yedges = y_km * 1e3

# process each datasets to get intensity map
inten_sets = []
radius_list = []
dataset_days = []
for file in filenames:
    # Read particle data
    with open(file, 'rb') as f:
        data = pickle.load(f)
        p_t = data['p_t']
        radius_dust = data['radii_dust']
        day = data['day']

    radius_list.append(radius_dust); dataset_days.append(day)

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

    # rotate particles with an angle on the projected plane
    # This step is to align simulated tail with observation in a line
    p_projected[:, 1], p_projected[:, 2] = rotate_coords(p_projected[:, 1], p_projected[:, 2], 7.2)
    
    """2D density histogram"""
    hist_data = p_projected[4:]

    # Extract x and y projected coordinates
    x_proj = hist_data[:, 1]
    y_proj = hist_data[:, 2]

    # Mask out particles outside of range
    mask = (
        (x_proj >= xedges[0]) & (x_proj <= xedges[-1]) &
        (y_proj >= yedges[0]) & (y_proj <= yedges[-1])
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
    # W = 1e6
    # weight = 10 ** ((np.log10(radius_dust)+1) * np.log10(W) / (-3))    
    weight = 1

    # scattering intensity (assume scattering phase angle constant for all particles)
    qext, qsca, qback, g = mie.efficiencies(m, 2*radius_dust, lambda0)
    p_func = 1
    px_inten = qsca * (np.pi * radius_dust**2) * px_den * weight * p_func
    px_inten = px_inten.T

    inten_sets.append(px_inten)

radius_list = np.array(radius_list)
dataset_days = np.array(dataset_days)
sim_stack = np.stack(inten_sets, axis=-1)            # shape: (ny, nx, n_sizes)

# set fitting input X and Y
fit_input_X = sim_stack.reshape(-1, len(filenames))  # shape: (ny*nx, n_sizes)
fit_input_Y = hst_data.flatten()                     # shape: (ny*nx)

# create a mask of HST meaningful signal (pixels of background signal 1e-20 are not used for fitting)
meaningful_signal_mask = np.log10(fit_input_Y + 1e-20) > -20.0

# get the meaningful pixels from the HST data and simulated data using the mask
fit_input_Y_filtered = fit_input_Y[meaningful_signal_mask]
fit_input_X_filtered = fit_input_X[meaningful_signal_mask, :]

# --- min/max values ---
sim_min = np.min(fit_input_X_filtered)
sim_max = np.max(fit_input_X_filtered)
hst_min = np.min(fit_input_Y_filtered)
hst_max = np.max(fit_input_Y_filtered)

# --- Formatted Print Statement ---
# The f-string formatting aligns the text and numbers into clean columns.
# :<20 means left-align in a 20-character space.
# :>15.2e means right-align in a 15-character space, formatted in scientific notation with 2 decimals.
print(f"{'':<20} {'min':>15} {'max':>15}")
print("-" * 55) # Optional: adds a separator line for clarity
print(f"{'Simulation':<20} {sim_min:>15.2e} {sim_max:>15.2e}")
print(f"{'HST Observation':<20} {hst_min:>15.2e} {hst_max:>15.2e}")

"""
Least Square Fitting
"""
from scipy.optimize import least_squares

def residuals(weights, I_models, I_obs):
    # weights can's be negative
    if np.any(weights < 0):
        return np.inf  # return a large error

    # calculate total intensity of the model
    I_fit = np.dot(I_models, weights)

    # calculate residuals in Logarithmic space 
    epsilon = 1e-20
    log_diff = np.log10(I_obs + epsilon) - np.log10(I_fit + epsilon)
    return log_diff

# guess an intial weights (e.g. all ones)
#initial_weights = np.ones(fit_input_X_filtered.shape[1])
initial_weights = np.array([1.48607575e-01, 2.03574515e-01, 3.44638220e-01, 4.84672084e-01,
                            5.36902768e-01, 5.75193493e-01, 6.08559074e-01, 5.99757568e-01,
                            5.61236310e-01, 4.91978292e-01, 3.86263984e-01, 2.45922593e-01,
                            7.13588707e-09])

# conduct the fitting
result = least_squares(
		residuals,
		initial_weights,
		args=(fit_input_X_filtered, fit_input_Y_filtered)
)

# final results
final_weights = result.x
#print("Fitted weights (particle number multipliers):\n", final_weights)

output_table = np.vstack([dataset_days, radius_list, final_weights]).T
np.savetxt(os.path.join(output_dir, 'fitting_output_83.77.txt'), output_table, fmt='%.8e', delimiter=',')
