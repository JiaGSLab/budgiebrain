import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load data: Using specific mm values for file names
from scipy import sparse
from scipy.stats import ttest_1samp
import os
file_path = "/data/output/"
mm_values = [0.3, 1, 3, 5, 6, 7, 8, 9, 10, 11]
scanpy_objects = [sc.read_h5ad(os.path.join(file_path, f"{mm}mm_Tel_RCTD.h5ad")) for mm in mm_values]

# Define gene list
geneset = ("CCN2","DACH2","DACH1", "SATB2", "SATB1", 
           "LYPD1", "SST", "PCP4", "MEIS2","MBP","ZNF385B",
           "SLC17A6","GAD2","NTS","PVALB","NR2F2","NRP2",
           'NRIP3', "CCK","SLC6A7","SLIT1","CACNA1H","GRIA4","BDNF","RARB"
          )
# 2. Find intersection of all genes across datasets
all_gene_sets = [set(obj.var_names) for obj in scanpy_objects]
common_genes = list(set.intersection(*all_gene_sets))

# Ensure the selected geneset is within the common genes
filtered_geneset = [gene for gene in geneset if gene in common_genes]
if not filtered_geneset:
    raise ValueError("None of the genes in geneset are present in all datasets.")
    

# 3. Calculate spatial scores and standard deviation for each slice
scores_per_slice = []
std_per_slice = []  # To store standard deviation for each slice
p_values_per_slice = []  # To store p-values for each gene in each slice
for mm, obj in zip(mm_values, scanpy_objects):
    obj = obj[:, common_genes]  # Align genes across datasets
    spatial_coords = obj.obsm['spatial']

    # Validate and process obj.X for compatibility
    if sparse.issparse(obj.X):
        expression_data = obj.X.todense()
    elif isinstance(obj.X, np.ndarray):
        expression_data = obj.X
    else:
        raise TypeError("Unsupported type for obj.X. Expected sparse or ndarray.")

    slice_df = pd.DataFrame(
        np.hstack([spatial_coords, expression_data]),
        columns=['x', 'y'] + common_genes
    )
    slice_df = slice_df[['x', 'y'] + filtered_geneset]  # Filter to selected genes

    # Convert gene expression to binary (0/1)
    threshold = 1  # Set threshold for binary conversion
    for gene in filtered_geneset:
        slice_df[gene] = (slice_df[gene] > threshold).astype(int)

    # Normalize coordinates for dorsal-ventral and lateral-medial directions
    slice_df['y_normalized'] = (slice_df['y'] - slice_df['y'].min()) / (slice_df['y'].max() - slice_df['y'].min())
    slice_df['x_normalized'] = (slice_df['x'] - slice_df['x'].min()) / (slice_df['x'].max() - slice_df['x'].min())

    # Calculate dorsal-ventral and lateral-medial scores and standard deviation
    dorsal_ventral_scores = []
    lateral_medial_scores = []
    dv_std_values = []
    lm_std_values = []
    p_values = []
    for gene in filtered_geneset:
        gene_spots = slice_df[slice_df[gene] > 0]
        if len(gene_spots) > 1:  # Ensure there are enough spots for calculation
            dv_std = np.std(gene_spots['y_normalized'])
            lm_std = np.std(gene_spots['x_normalized'])
            dv_score = np.mean(gene_spots['y_normalized'])
            lm_score = np.mean(gene_spots['x_normalized'])
            t_stat, p_val = ttest_1samp(gene_spots['y_normalized'], popmean=0.5)
        else:
            dv_std = lm_std = dv_score = lm_score = p_val = np.nan

        dorsal_ventral_scores.append(dv_score)
        lateral_medial_scores.append(lm_score)
        dv_std_values.append(dv_std)
        lm_std_values.append(lm_std)
        p_values.append(p_val)

    # Combine into a dataframe for the slice
    scores_per_slice.append(pd.DataFrame({
        'Gene': filtered_geneset,
        'Dorsal_Ventral_Score': dorsal_ventral_scores,
        'Lateral_Medial_Score': lateral_medial_scores,
        'DV_Std': dv_std_values,
        'LM_Std': lm_std_values,
        'P_Value': p_values,
        'Slice': [mm] * len(filtered_geneset)
    }))

# Combine all slices into one dataframe
scores_df = pd.concat(scores_per_slice, ignore_index=True)

# Ensure no duplicate entries by resetting index on pivoting
scores_df = scores_df.drop_duplicates()

# 4. Plot dorsal-ventral heatmap 
heatmap_data_dv = scores_df.pivot(index='Gene', columns='Slice', values='Dorsal_Ventral_Score')
plt.figure(figsize=(12, 8))
sns.heatmap(
    heatmap_data_dv,
    cmap='coolwarm',
    annot=True,
    cbar_kws={'label': 'Dorsal-Ventral Score'}
)
plt.title('Dorsal-Ventral Scores per Slice')
plt.xlabel('Slice (mm)')
plt.ylabel('Gene')
plt.show()
heatmap_data_dv.to_csv('/data/work/spatial_distribution/heatmap_data_dv.csv', index=False)


# 5. Plot lateral-medial heatmap
heatmap_data_lm = scores_df.pivot(index='Gene', columns='Slice', values='Lateral_Medial_Score')
plt.figure(figsize=(12, 8))
sns.heatmap(
    heatmap_data_lm,
    cmap='coolwarm',
    annot=True,
    cbar_kws={'label': 'Lateral-Medial Score'}
)
plt.title('Lateral-Medial Scores per Slice')
plt.xlabel('Slice (mm)')
plt.ylabel('Gene')
plt.show()
heatmap_data_lm.to_csv('/data/work/spatial_distribution/heatmap_data_lm.csv', index=True)