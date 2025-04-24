import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
adata = sc.read_h5ad('/data/work/figures/annotation/5.h5ad')
rctd_results = pd.read_csv('/data/work/RCTD/5mm/5mm_norm_weights.csv', index_col=0)
rctd_results.index = rctd_results.index.astype(str)
common_indices = adata.obs.index.intersection(rctd_results.index)
rctd_results_filtered = rctd_results.loc[common_indices]
adata_filtered = adata[common_indices, :]
adata_filtered.obs = adata_filtered.obs.join(rctd_results_filtered, how='left')
adata = adata_filtered
#plot example
sc.pl.spatial(
    adata,
    color='MGE_SST_PVALB_Inh',  # su
     cmap='magma',   # 使用更适合分类数据的颜色调色板
    spot_size=150,        # 调整点的大小
    show=True
)

# 显示图像
plt.show()