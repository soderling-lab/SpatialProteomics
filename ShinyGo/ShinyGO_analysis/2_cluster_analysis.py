
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics.pairwise import euclidean_distances
import os
# Load the data
currdir = os.getcwd()
parent = os.path.dirname(currdir)
print("PARENT: ", parent)
GENE_NAME = "LRRK2"
data = pd.read_csv(f'{parent}/ShinyGO_data/raw_heatmaps/LRRK2_heatmap_top_pathways_BiologicalComponent.csv')
PATHWAY = "BIOLOGICALcmpnt"
# Replace -1.000000e-09 with NaN
data = data.replace(-1.000000e-09, np.nan)

# Extract module columns
modules = data.iloc[:, 1:]

# Replace NaN values with a large number before computing distances
# We use a large number to simulate 'infinite' distance in missing data points
modules_filled = modules.fillna(1e10)

# Compute the distance matrix
dist_matrix = euclidean_distances(modules_filled, modules_filled)

# Perform hierarchical clustering
link = linkage(dist_matrix, 'ward')

# Get the order of rows according to the hierarchy
order = leaves_list(link)

# Reorder the data according to the clustering result
data_reordered = data.iloc[order, :]

# Save the reordered data to a new CSV file
data_reordered.to_csv(f'{parent}/ShinyGO_data/rendered_heatmaps/{GENE_NAME}_heatmap_LRRK2_pathways_reordered_{PATHWAY}.csv', index=False)

print("saved")
