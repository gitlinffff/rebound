"""
shift particle coordinates in data_high.txt 
"""

import numpy as np
import os

def shift_y_coordinate(input_file, output_file, shift_distance):
    """
    Reads particle data, shifts the y-coordinate, and saves to a new file.

    Args:
        input_file (str): Path to the source data file.
        output_file (str): Path to save the modified data file.
        shift_distance (float): The distance to shift the y-coordinates (in cm).
    """
    print(f"\nReading data from '{input_file}'...")

    # Load the numerical data, skipping the header row
    data = np.loadtxt(input_file)

    print(f"Shifting y-coordinates by {shift_distance} cm...")
    # The y-coordinate column
    data[:, 2] += shift_distance

    print(f"Saving new data to '{output_file}'...")
    # Save the modified array to the output file, preserving the format and header
    np.savetxt(output_file, data, fmt='%.18e', delimiter=' ')
    print("Done.")


if __name__ == '__main__':
    # Define file names and the shift distance
    input_filename = "/nuke/linfel/Ejecta/data_high.txt"
    output_filename = "/nuke/linfel/Ejecta/data_high_shifted.txt"
    y_shift = -5500.0  # The amount to move particles along the y-axis (in cm)

    shift_y_coordinate(input_filename, output_filename, y_shift)

    # 3. Verify the result
    print("\n--- Verification ---")
    original_data = np.loadtxt(input_filename)
    shifted_data = np.loadtxt(output_filename)

    print("Original Y-coords | Shifted Y-coords")
    for i in range(5): # Print the first 5 rows for comparison
        print(f"{original_data[i, 2]:>17.2f} | {shifted_data[i, 2]:>17.2f}")
