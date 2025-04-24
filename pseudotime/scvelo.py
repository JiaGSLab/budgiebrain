#scvelo 0.2.5 (python 3.8.15)
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
adata_combined = sc.read_h5ad('/data/work/Tel_dev/scVelo/data/scvelo_merged.h5ad')
cell_mask = adata_combined.obs['dev_celltype'].isin(['Meso_GAS2_Ipc',  'Meso_KIAA1217_Ex',
                                                    'Meso_RARB_Ex','Pallium_Ipc','Meso_SOX6_Immature_Ex',
                                                     'Progenitor_Ex', 'Meso_SOX6_Ex'])

Meso_data = adata_combined[cell_mask]
scv.pp.filter_and_normalize(adata_combined, min_shared_counts=20, n_top_genes=3000)
scv.pp.moments(adata_combined, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_combined)
scv.tl.velocity_graph(adata_combined)
scv.pl.velocity_embedding_stream(adata_combined, basis='umap_harmony',color='dev_celltype')