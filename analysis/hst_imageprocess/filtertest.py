import os
from coordinates import day_hstfile_mapping
from tail_intensity_2 import load_and_preprocess_fits
from process_hst import select_area_and_get_indices, plot_selected_region


if __name__ == "__main__":
# Define file paths and coordinate key
	day_code = 'day_64.44' # Key for coordinate data
	FILEPATH_1 = os.path.join("/home/linfel/linfel_data/hst_raw_JianyangLi/", day_hstfile_mapping[day_code])
	output_dir = f"/home/linfel/linfel_data/longterm_anal/{day_code}_380-450"

# Load and process Image 1 (This sets the primary coordinate system)
	log10_hst, extent, pixel_km, nx, ny = load_and_preprocess_fits(FILEPATH_1, day_code)
	polygon_mask_1d, vertices = select_area_and_get_indices(log10_hst, extent, nx, ny,
	                                          mask_filepath=os.path.join(output_dir, "hst_region_1dmask.npy"))
	#plot_hst_image_with_selection(log10_hst, extent, vertices,  output_dir, day_code)
	plot_selected_region(log10_hst, polygon_mask_1d, extent, output_dir, day_code)

