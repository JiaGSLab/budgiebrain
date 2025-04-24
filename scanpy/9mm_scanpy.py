import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/9mm.h5ad')
# annotation
simplified_cluster_annotations = {
 1: "ZIC1+ ZIC2+ Nucleus of Thalamus",
    2: "Thalamus",
    3: "hyper-nidopallium",
    4: "Mesopallium",
    5: "SGF",
    6: "NIVL",
    7: "Field L",
    8: "hyper-nidopallium",
    9: "CHGB+ CALB2+ Nucleus of Thalamus",
    10: "NIVL",
    11: "Arcopallium",
    12: "Imc/Ipc",
    13: "hyper-nidopallium",
    14: "hyper-nidopallium",
    15: "White matter",
    16: "Thalamus",
    17: "hyper-nidopallium",
    18: "Striatum",
    19: "Mesopallium",
    20: "Vasculature",
    21: "Hippocampus",
    22: "SGP",
    23: "Vasculature",
    24: "Mesopallium",
    25: "Pons",
    26: "ZIC1+ ZIC2+ Nucleus of Thalamus",
    27: "SOp",
    28: "Vasculature",
    29: "NA",
    30: "ZIC1+ ZIC2+ Nucleus of Thalamus",
    31: "Nidopallium Intermediale",
    32: "SGC",
    33: "Arcopallium",
    34: "Vasculature",
    35: "Nidopallium Intermediale",
    36: "Nidopallium Intermediale",
    37: "NLc"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)

# split "hyper-nidopallium"
spatial_coords = adata.obsm['spatial']
clusters = adata.obs['simplified_annotated_spatial_leiden']
new_annotations = clusters.copy()
target_cluster = "hyper-nidopallium"
existing_categories = clusters.cat.categories.tolist()
additional_categories = ["Hyperpallium", "Nidopallium Intermediale"]
categories_to_add = [cat for cat in additional_categories if cat not in existing_categories]
if categories_to_add:
    new_annotations = new_annotations.cat.add_categories(categories_to_add)
hyper_nidopallium_mask = (clusters == target_cluster)
hyperpallium_mask = hyper_nidopallium_mask & (spatial_coords[:, 0] > 15000) & (spatial_coords[:, 1] > 16000)
new_annotations[hyperpallium_mask] = "Hyperpallium"
new_annotations[hyper_nidopallium_mask & ~hyperpallium_mask] = "Nidopallium Intermediale"
adata.obs['simplified_annotated_spatial_leiden'] = new_annotations
print(adata.obs['simplified_annotated_spatial_leiden'].value_counts())

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
    'Hippocampus':'Hp',
    'SGP':'SGP',
    'SOp':'SOp',
    'SGF':'SGF',
    'SGC':'SGC',
    'Field L':'Field L',
    'Imc/Ipc':'Imc/Ipc',
    'Pons':'BS',
    'NA':'NA'
}

current_annotations = adata.obs['simplified_annotated_spatial_leiden']
updated_annotations = current_annotations.map(simplified_cluster_annotations)
adata.obs['simplified_annotated_spatial_leiden'] = updated_annotations

#filter NA
valid_cells = adata.obs['simplified_annotated_spatial_leiden'] != 'NA'
adata_filtered = adata[valid_cells].copy()
print("\n filtered AnnData:")
print(adata_filtered)