import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/8mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Str",
    2: "NI",
    3: "Meso",
    4: "NI",
    5: "CHGB+ CALB2+ nTh",
    6: "Th",
    7: "WM",
    8: "WM",
    9: "NIVL",
    10: "SGF",
    11: "ZIC1+ ZIC2+ nTh",
    12: "Th",
    13: "NI",
    14: "SGF",
    15: "Vasc",
    16: "SOp",
    17: "HSPB7+ nTh",
    18: "NI",
    19: "NIVL",
    20: "Arco",
    21: "Hyper",
    22: "WM",
    23: "Arco",
    24: "Ento",
    25: "SGC"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)