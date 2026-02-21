import os
import numpy as np
import matplotlib.pyplot as plt
from hst_fitting import load_array_from_h5
from coordinates import position_dict

def get_hubble_pixel_km(day_code):
	"""
	get the Hubble pixel size at the target's location.
	"""
# Calculate target distance
	r_Didy_sys_bary = position_dict['Didy_sys_bary'][day_code] # Sun Body Center frame
	r_Hubble = position_dict['Hubble'][day_code]
	target_distance = np.linalg.norm(r_Didy_sys_bary - r_Hubble)  # km

# Hubble pixel size
	pixel_arcsec = 0.04
	range_km = target_distance
	pixel_fov = np.deg2rad(pixel_arcsec / 3600)  # pixel field of view in rad
	pixel_km = 2 * range_km * np.tan(pixel_fov/2)

	return pixel_km

def plot_fitted_image(I_fit, pixel_km, vmin, output_dir):
  # axes scale
  ny, nx = I_fit.shape
  x_km = np.arange(nx+1) * pixel_km
  y_km = np.arange(ny+1) * pixel_km

  # Create meshgrid for bin edges
  X, Y = np.meshgrid(x_km, y_km)   # km

  # Plot pcolormesh
  plt.figure(figsize=(8, 8))

  log10_fit = np.log10(I_fit + 1e-20)
  pc = plt.pcolormesh(X, Y, log10_fit, cmap='cividis', shading='auto', vmin=vmin, vmax=np.nanmax(log10_fit))

  plt.xlabel('Projected X [km]')
  plt.ylabel('Projected Y [km]')
  plt.title('Intensity Fitting on View Plane')
  plt.grid(True, linestyle='--', linewidth=0.1, color='red', alpha=0.7)
  plt.gca().set_aspect('equal', adjustable='box')

  # Colorbar
  cbar = plt.colorbar(pc, orientation='horizontal', pad=0.1, shrink=0.5, aspect=30)
  cbar.set_label(r'$\log_{10}$(Brightness) [$W m_{-2} um_{-1} sr_{-1}$]')

  plt.tight_layout()
  output_name = os.path.join(output_dir, f"fitted_image_vmin{vmin}.png")
  plt.savefig(output_name, dpi=300, bbox_inches='tight', pad_inches=0.1)
  #plt.show()
  plt.close()


def main():
	day_code = "day_14.91"
	wdir = "/home/linfel/linfel_data/shortterm_anal/day_14.91_fabio_radial_discrete/day_14.91_270-410_step8"
	pixel_km = get_hubble_pixel_km(day_code)
	I_fit = load_array_from_h5(os.path.join(wdir, 'I_fit.h5'), 'intensity')

	vmin = [-7,-6,-5,-4,-3]
	vmin = [-8]
	for value in vmin:
		plot_fitted_image(I_fit, pixel_km, value, wdir)


if __name__ == "__main__":
	main()
