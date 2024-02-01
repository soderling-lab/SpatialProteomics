
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics.pairwise import euclidean_distances
import os
import nbimporter
from module_pathways import ensure_dirs_exists

# Load the data
currdir = os.getcwd()
parent = os.path.dirname(currdir)
gparent = os.path.dirname(parent)
GENE_NAME = "LRRK2"
PATHWAY = "BiologicalComponent"

data = pd.read_csv(f'{parent}/example_data/ShinyGO_data/raw_heatmaps/{GENE_NAME}/{GENE_NAME}_heatmap_top_pathways_BiologicalComponent.csv')
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

SAVE_FILE = f'{parent}/example_data/ShinyGO_data/rendered_heatmaps/{GENE_NAME}/{PATHWAY}_heatmap_pathways_reordered.csv'
ensure_dirs_exists(SAVE_FILE)
data_reordered.to_csv(SAVE_FILE, index=False)

print("saved")
