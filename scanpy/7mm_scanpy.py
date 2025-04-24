import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/7mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Striatum",
    2: "White matter",
    3: "Mesopallium",
    4: "hyper-nidopallium",
    5: "Nidopallium Intermediale",
    6: "White matter",
    7: "NIVL",
    8: "CHGB+ CALB2+ Nucleus of Thalamus",
    9: "Mesopallium",
    10: "Arcopallium",
    11: "ZIC1+ ZIC2+ Nucleus of Thalamus",
    12: "NLc",
    13: "Arcopallium",
    14: "Mesopallium",
    15: "Hyperpallium",
    16: "hyper-nidopallium",
    17: "hyper-nidopallium",
    18: "Vasculature",
    19: "Nidopallium Intermediale",
    20: "CHGB+ CALB2+ Nucleus of Thalamus",
    21: "Thalamus",
    22: "hyper-nidopallium",
    23: "Mesopallium",
    24: "Intercalated Nidopallium",
    25: "hyper-nidopallium",
    26: "White matter",
    27: "Thalamus",
    28: "Thalamus",
    29: "Vasculature",
    30: "NIVL",
    31: "Hyperpallium",
    32: "GBX2+ HSPB7+ Nucleus of Thalamus",
    33: "Nidopallium Intermediale"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)

# split "hyper-nidopallium"
spatial_coords = adata.obsm['spatial']
clusters = adata.obs['simplified_annotated_spatial_leiden']
new_annotations = clusters.copy()
target_cluster = "hyper-nidopallium"
new_annotations[(clusters == target_cluster) & (spatial_coords[:, 1] > 6200)] = "Nidopallium Intermediale"
new_annotations[(clusters == target_cluster) & (spatial_coords[:, 1] <= 6200)] = "Hyperpallium"
adata.obs['simplified_annotated_spatial_leiden'] = new_annotations

# update_annotation
simplified_cluster_annotations = {
    'Striatum': 'Str',
    'Nidopallium Intermediale': 'NI',
    'Mesopallium': 'Meso',
    'CHGB+ CALB2+ Nucleus of Thalamus': 'CHGB+ CALB2+ nTh',
    'Thalamus': 'Th',
    'White matter': 'WM',
    'NIVL': 'NIVL',
    'Vasculature': 'Vasc',
    'Arcopallium': 'Arco',
    'Hyperpallium': 'Hyper',
    'Intercalated Nidopallium': 'Ento',
    'ZIC1+ ZIC2+ Nucleus of Thalamus': 'ZIC1+ ZIC2+ nTh',
    'GBX2+ HSPB7+ Nucleus of Thalamus': 'HSPB7+ nTh',
}

current_annotations = adata.obs['simplified_annotated_spatial_leiden']
updated_annotations = current_annotations.map(simplified_cluster_annotations)
adata.obs['simplified_annotated_spatial_leiden'] = updated_annotations