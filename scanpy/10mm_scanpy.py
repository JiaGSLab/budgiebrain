import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/10mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "Th",
    2: "NI",
    3: "Meso",
    4: "CGL",
    5: "NC",
    6: "Vasc",
    7: "NI",
    8: "CML",
    9: "SGF",
    10: "CPL",
    11: "SGF",
    12: "Field L",
    13: "SGP",
    14: "CHGB+ CALB2+ nTh",
    15: "CML",
    16: "NI",
    17: "Hp",
    18: "NI",
    19: "WM",
    20: "Th",
    21: "SGC",
    22: "SOp",
    23: "Imc/Ipc",
    24: "Vasc",
    25: "CML",
    26: "CML",
    27: "Hp",
    28: "WM",
    29: "Th",
    30: "Vasc"
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)