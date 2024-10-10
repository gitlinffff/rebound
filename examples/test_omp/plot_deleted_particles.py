import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the data from a CSV file (without header)
# Adjust the file name/path accordingly
data = pd.read_csv('deleted_particles.csv', header=None)

# Manually assign column names: x, y, z, category
data.columns = ['x', 'y', 'z', 'category']

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Define a colormap for different categories
categories = data['category'].unique()
colors = plt.cm.get_cmap('Set1', len(categories))

# Plot particles with different colors for different categories
for i, category in enumerate(categories):
    if category == 3: continue
    category_data = data[data['category'] == category]
    ax.scatter(category_data['x'], category_data['y'], category_data['z'],
               color=colors(i), label=f'Category {category}', s=50)
# Add labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Scatter Plot of Particles')

# Add legend
ax.legend()

# Show the plot
plt.savefig("p1.png",dpi=150)

