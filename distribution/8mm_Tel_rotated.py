#8mm
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
adata = sc.read_h5ad('/data/work/figures/annotation/8mm.h5ad')
rctd_results = pd.read_csv('/data/work/RCTD/8mm/8mm_norm_weights.csv', index_col=0)
rctd_results.index = rctd_results.index.astype(str)
common_indices = adata.obs.index.intersection(rctd_results.index)
rctd_results_filtered = rctd_results.loc[common_indices]
adata_filtered = adata[common_indices, :]
adata_filtered.obs = adata_filtered.obs.join(rctd_results_filtered, how='left')
adata = adata_filtered
spatial_coords = adata.obsm['spatial']

plt.figure(figsize=(6, 8))
plt.scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='blue', label='Original', alpha=0.5)
plt.title("Original Spatial Coordinates")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.show()

#rotation
theta =90  
theta_rad = np.deg2rad(theta)  
rotation_matrix = np.array([[np.cos(theta_rad), -np.sin(theta_rad)],
                            [np.sin(theta_rad), np.cos(theta_rad)]])


spatial_coords = adata.obsm['spatial']
rotated_coords = spatial_coords @ rotation_matrix.T  
adata.obsm['spatial_rotated'] = rotated_coords
adata.obsm['spatial'] = adata.obsm['spatial_rotated']
spatial_coords = adata.obsm['spatial']

# flipped
flipped_coords = spatial_coords.copy()  
flipped_coords[:, 0] = -flipped_coords[:, 0]
adata.obsm['spatial_flipped'] = flipped_coords
adata.obsm['spatial'] = adata.obsm['spatial_flipped']

#fliter
spatial_coords = adata.obsm['spatial']
mask = spatial_coords[:, 1] >= 5000 
adata = adata[mask, :]  

#save
adata.write("/data/output/8mm_Tel_RCTD.h5ad")