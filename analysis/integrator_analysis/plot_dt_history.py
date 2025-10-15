import numpy as np
import os
import matplotlib.pyplot as plt

# List of paths to your dt_history.csv files
file_paths = [
	'/home/linfel/linfel_scratch/rebound_exp/test_integrator/bs_um_004/dt_history.csv',
	'/home/linfel/linfel_scratch/rebound_exp/test_integrator/ias_mm_005/dt_history.csv',
	'/home/linfel/linfel_scratch/rebound_exp/test_integrator/ias_mm_006/dt_history.csv',
	'/home/linfel/linfel_scratch/rebound_exp/test_integrator/bs_mm_008/dt_history.csv',
	#'/home/linfel/linfel_scratch/rebound_exp/test_integrator/ias_mm_005/dt_history.csv',
]

# output path
output_dir = '/home/linfel/linfel_scratch/rebound_exp/test_integrator/postprocess/'
os.makedirs(output_dir, exist_ok=True)

plt.figure(figsize=(10, 6))

for i, file_path in enumerate(file_paths):
	label = file_path.split('/')[-2]
	data = np.genfromtxt(file_path, skip_header=1)
	plt.plot(data[:,1], data[:,3], lw=0.5, label=label)

plt.yscale('log')
plt.xlabel('t')
plt.ylabel('dt_minimum')
plt.title('dt Time Series for Multiple Runs')
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_dir, 'dt_history_multiple_runs.png'), dpi=150, bbox_inches='tight', pad_inches=0.1)
plt.show()
