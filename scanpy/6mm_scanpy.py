import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/6mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Striatum",
    2: "Striatum",
    3: "Striatum",
    4: "NIVL",
    5: "Mesopallium",
    6: "hyper-nidopallium",
    7: "hyper-nidopallium",
    8: "Mesopallium",
    9: "Striatum",
    10: "hyper-nidopallium",
    11: "Nidopallium Intermediale",
    12: "Striatum",
    13: "Intercalated Nidopallium",
    14: "Striatum",
    15: "hyper-nidopallium",
    16: "RGS12+ MStm",
    17: "Vasculature",
    18: "White matter",
    19: "Hyperpallium",
    20: "Mesopallium",
    21: "Striatum",
    22: "Striatum",
    23: "Striatum",
    24: "Vasculature",
    25: "Vasculature",
    26: "Vasculature",
    27: "Mesopallium",
    28: "Nidopallium Intermediale",
    29: "NIVL",
    30: "Striatum",
    31: "Vasculature",
    32: "Striatum"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)

# split "hyper-nidopallium"
spatial_coords = adata.obsm['spatial']
clusters = adata.obs['simplified_annotated_spatial_leiden']
new_annotations = clusters.copy()
target_cluster = "hyper-nidopallium"
print("Initial hyper-nidopallium count:", np.sum(clusters == target_cluster))
hyperpallium_mask = (clusters == target_cluster) & (spatial_coords[:, 0] < 9100) & (spatial_coords[:, 1] < 19000)
nidopallium_intermediale_mask = (clusters == target_cluster) & ((spatial_coords[:, 0] >= 9100) | (spatial_coords[:, 1] >= 19000))
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
    'NIVL': 'NIVL',
    'Hyperpallium': 'Hyper',
    'Intercalated Nidopallium': 'Ento',
    "RGS12+ MStm": 'RGS12+ MStm'
}

current_annotations = adata.obs['simplified_annotated_spatial_leiden']
updated_annotations = current_annotations.map(simplified_cluster_annotations)
adata.obs['simplified_annotated_spatial_leiden'] = updated_annotations