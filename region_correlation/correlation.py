import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# set path
data_path = '/data/work/figures/annotation/'

# loading
h5ad_files = [f for f in os.listdir(data_path) if f.endswith('.h5ad')]
adatas = []
for f in h5ad_files:
    adata = sc.read(os.path.join(data_path, f))
    adata.obs['batch'] = os.path.splitext(f)[0]  
    adatas.append(adata)

# merge
adata_combined = adatas[0].concatenate(adatas[1:], batch_key='batch', batch_categories=[os.path.splitext(f)[0] for f in h5ad_files])
print(adata_combined)
output_file = os.path.join(data_path, 'merged_data.h5ad')
adata_combined.write(output_file)

#loading merge data
adata = sc.read_h5ad('/data/work/figures/annotation/merged_data.h5ad')

# mapping
annotation_map = {
    'Hyper': ['Hyper'],
    'Meso': ['Meso'],
    'Nido': ['NF', 'NI', 'NIVL', 'Field L', 'NCM', 'NCC', 'NCL', 'NC'],
    'Ento': ['Ento'],
    'Arco': ['Arco'],
    'Hp': ['Hp'],
    'Str': ['Str', 'RGS12+ MStm'],
    'Th': ['CHGB+ CALB2+ nTh', 'Th', 'ZIC1+ ZIC2+ nTh'],
    'Cere': ['CGL', 'CPL'],
    'TeO': ['SGF', 'SGC', 'SGP', 'Imc/Ipc','TeO'],
    'Others': ['SOp', 'WM', 'BS', 'Vasc', 'CML', 'HSPB7+ nTh']
}

reversed_annotation_map = {item: key for key, values in annotation_map.items() for item in values}
adata.obs['annotation'] = adata.obs['simplified_annotated_spatial_leiden'].map(reversed_annotation_map)
print(adata.obs[['simplified_annotated_spatial_leiden', 'annotation']])

# Calculate the mean expression per cluster
cluster_means = adata.to_df().groupby(adata.obs['simplified_annotated_spatial_leiden']).mean()

# Compute the Pearson correlation matrix
correlation_matrix = cluster_means.T.corr(method='pearson')

# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Cluster Pearson Correlation Heatmap')
plt.show()