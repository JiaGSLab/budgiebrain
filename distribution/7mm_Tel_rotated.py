#0.3mm
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
adata = sc.read_h5ad('/data/work/figures/annotation/7mm.h5ad')
spatial_coords = adata.obsm['spatial']

plt.figure(figsize=(6, 8))
plt.scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='blue', label='Original', alpha=0.5)
plt.title("Original Spatial Coordinates")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.show()


spatial_coords = adata.obsm['spatial']

# flipped
flipped_coords = spatial_coords.copy()  
flipped_coords[:, 1] = -flipped_coords[:, 1] 
adata.obsm['spatial_flipped'] = flipped_coords
adata.obsm['spatial'] = adata.obsm['spatial_flipped']

#fliter
spatial_coords = adata.obsm['spatial']
mask = spatial_coords[:, 0] >= 5000
adata = adata[mask, :]  

#save
adata.write("/data/output/7mm_Tel_RCTD.h5ad")