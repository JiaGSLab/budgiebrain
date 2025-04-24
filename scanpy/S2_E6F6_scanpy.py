import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/S2_E6F6.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Str",
    2: "NI",
    3: "Meso",
    4: "Th",
    5: "CGL",
    6: "Hyper",
    7: "CML",
    8: "Str",
    9: "Hyper",
    10: "NF",
    11: "Vasc",
    12: "Meso",
    13: "Ento",
    14: "Th",
    15: "WM",
    16: "BS",
    17: "CPL",
    18: "NC",
    19: "NF",
    20: "Field L",
    21: "WM",
    22: "HSPB7+ nTh",
    23: "Hp",
    24: "Vasc",
    25: "Str",
    26: "NF",
    27: "CHGB+ CALB2+ nTh"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)