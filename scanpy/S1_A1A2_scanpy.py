import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/S.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Th",
    2: "NI",
    3: "Str",
    4: "Hyper",
    5: "Meso",
    6: "WM",
    7: "NF",
    8: "NI",
    9: "Vasc",
    10: "CGL",
    11: "Str",
    12: "Hyper",
    13: "Ento",
    14: "BS",
    15: "CPL",
    16: "NC",
    17: "Hp",
    18: "Meso",
    19: "CML",
    20: "CHGB+ CALB2+ nTh",
    21: "WM",
    22: "BS",
    23: "Th",
    24: "HSPB7+ nTh",
    25: "Str",
    26: "Field L",
    27: "TeO",
    28: "Ento"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)