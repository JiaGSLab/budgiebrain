import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# loading
adata = sc.read_h5ad('/data/work/figures/annotation/1mm.h5ad')
# annotation
simplified_cluster_annotations = {
    1: "NF",
    2: "NF",
    3: "Meso",
    4: "Meso",
    5: "Str",
    6: "NF",
    7: "NF",
    8: "Vasc",
    9: "NF",
    10: "Hyper",
    11: "NF"
 
}


# add_info
adata.obs['simplified_annotated_spatial_leiden'] = adata.obs['spatial_leiden'].astype(int).map(simplified_cluster_annotations)