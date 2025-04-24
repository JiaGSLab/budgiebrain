import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/5mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Mesopallium",
    2: "Striatum",
    3: "RGS12+ MStm",
    4: "hyper-nidopallium",
    5: "hyper-nidopallium",
    6: "Striatum",
    7: "Striatum",
    8: "Hyperpallium",
    9: "Striatum",
    10: "Mesopallium",
    11: "Striatum",
    12: "Vasculature",
    13: "Striatum",
    14: "hyper-nidopallium",
    15: "Nidopallium Intermediale",
    16: "hyper-nidopallium",
    17: "Intercalated Nidopallium",
    18: "hyper-nidopallium",
    19: "White matter",
    20: "Hyperpallium",
    21: "Striatum",
    22: "Striatum",
    23: "Vasculature",
    24: "Nidopallium Intermediale"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)

# split "hyper-nidopallium"
spatial_coords = adata.obsm['spatial']
clusters = adata.obs['simplified_annotated_spatial_leiden']
new_annotations = clusters.copy()
target_cluster = "hyper-nidopallium"
print("Initial hyper-nidopallium count:", np.sum(clusters == target_cluster))
hyperpallium_mask = (clusters == target_cluster) & (spatial_coords[:, 1] <= 10500) & (spatial_coords[:, 0] < 13500)
nidopallium_intermediale_mask = (clusters == target_cluster) & ((spatial_coords[:, 1] > 10500) | (spatial_coords[:, 0] >= 13500))
print("Hyperpallium mask count:", np.sum(hyperpallium_mask))
print("Nidopallium Intermediale mask count:", np.sum(nidopallium_intermediale_mask))
new_annotations[hyperpallium_mask] = "Hyperpallium"
new_annotations[nidopallium_intermediale_mask] = "Nidopallium Intermediale"
adata.obs['simplified_annotated_spatial_leiden'] = new_annotations
updated_counts = adata.obs['simplified_annotated_spatial_leiden'].value_counts()
print(updated_counts)

# update_annotation
simplified_cluster_annotations = {
    'Striatum': 'Str',
    'Nidopallium Intermediale': 'NI',
    'Mesopallium': 'Meso',
    'White matter': 'WM',
    'Vasculature': 'Vasc',
    'Hyperpallium': 'Hyper',
    'Intercalated Nidopallium': 'Ento',
    "RGS12+ MStm": 'RGS12+ MStm'
}

current_annotations = adata.obs['simplified_annotated_spatial_leiden']
updated_annotations = current_annotations.map(simplified_cluster_annotations)
adata.obs['simplified_annotated_spatial_leiden'] = updated_annotations