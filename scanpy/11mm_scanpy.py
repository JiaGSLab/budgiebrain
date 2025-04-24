import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/11mm.h5ad')
# annotation
simplified_cluster_annotations = {
 1: "WM",
    2: "NCC",#NCC
    3: "WM",
    4: "CML",
    5: "Vasc",
    6: "NCM",#NCM
    7: "BS",
    8: "SGF",
    9: "CGL",
    10: "CML",
    11: "BS",
    12: "BS",
    13: "NCL",#NCL
    14: "CML",
    15: "BS",
    16: "Th",
    17: "CML",
    18: "SGF",
    19: "SGP",
    20: "Vasc",
    21: "SGC",
    22: "CML",
    23: "CPL",
    24: "Imc",
    25: "Hp",
    26: "CML",
    27: "Th",
    28: "Vasc",
    29: "WM",
    30: "NCC",#NCC
    31: "BS",
    32: "CML",
    33: "WM"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)