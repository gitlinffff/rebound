import os
import pandas as pd
import matplotlib.pyplot as plt


def plot_w_rl(wdir):
	"""
	Plot fitted weights against dust sizes, seperated by different radial distance bins.
	"""
	
	filepath = os.path.join(wdir, "w_r.csv")
	if not os.path.exists(filepath):
				print(f"Error: {filepath} not found.")
				return

	# Load data
	df = pd.read_csv(filepath)

	# Plotting script logic
	plt.figure(figsize=(10, 6))

	# Grouping by 'radial_group'
	for name, group in df.groupby('radial_group'):
			# Sorting by radius ensures the line plots correctly
			group = group.sort_values('radius')
			
			# We plot radius vs weights with errors
			plt.errorbar(group['radius'], group['weights'], yerr=group['errors'], 
									 label=f'Radial Bin {int(name)}', fmt='-o', 
									 markersize=4, capsize=3, elinewidth=1, alpha=0.8)

	plt.xscale('log') # Useful if radius spans multiple orders
	plt.yscale('log') # Often weights/mass follow power laws
	plt.xlabel('Particle Radius [m]')
	plt.ylabel('Fitted Weight (Mass Scaling)')
	plt.title('Fitted Weights vs Grain Size for different Radial Bins')
	plt.grid(True, which="both", ls="-", alpha=0.5)
	plt.legend(title='Distance from Dimorphos', bbox_to_anchor=(1.05, 1), loc='upper left')
	plt.tight_layout()

	plt.savefig(os.path.join(wdir, 'weights_distribution_plot.png'))
	print("Plot saved as weights_distribution_plot.png")
	plt.close()

if __name__ == "__main__":
	for step in [8, 4, 2, 1]:
		for stop in [410, 450]:	
			wdir = f"/home/linfel/linfel_data/shortterm_anal/day_14.91_fabio_radial_discrete/day_14.91_270-{stop}_step{step}"
			plot_w_rl(wdir)
