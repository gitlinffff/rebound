import numpy as np
import matplotlib.pyplot as plt

data_8377 = np.genfromtxt('/home/linfel/linfel_data/ejecta_exp_datahigh/postprocess/weight_fit_result/fitting_output_83.77.txt', delimiter=',')
data_13129 = np.genfromtxt('/home/linfel/linfel_data/ejecta_exp_datahigh/postprocess/weight_fit_result/fitting_output_131.29.txt', delimiter=',')


plt.loglog(data_8377[:,1], data_8377[:,2], label=r"$T_0$+83.77")
plt.loglog(data_13129[:,1], data_13129[:,2], label=r"$T_0$+131.29")

# Labels and title
plt.xlabel("Radius")
plt.ylabel("Weights")
plt.title("Ejecta Weight Distribution")

# Legend and grid
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.7)

# Save high-quality PNG
plt.tight_layout()
plt.savefig("/home/linfel/linfel_data/ejecta_exp_datahigh/postprocess/weight_distribution.png", dpi=150)
