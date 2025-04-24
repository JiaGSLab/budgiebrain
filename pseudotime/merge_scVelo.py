import scanpy as sc
import os
data_path = '/data/work/Tel_dev/scVelo/data/'
h5ad_files = [f for f in os.listdir(data_path) if f.endswith('.h5ad')]
adatas = []
for f in h5ad_files:
    adata = sc.read(os.path.join(data_path, f))
    adata.obs['Sample_ID'] = os.path.splitext(f)[0]  
    adatas.append(adata)

# merge
adata_combined = adatas[0].concatenate(adatas[1:], batch_key='Sample_ID', batch_categories=[os.path.splitext(f)[0] for f in h5ad_files])

print(adata_combined)
adata_total =  sc.read_h5ad("/data/work/Tel_dev/scVelo/dev_Tel_annotated.h5ad")
adata_total
common_gene = list(adata_total.var_names.intersection(adata_combined.var_names))
adata_total = adata_total[:, common_gene]
adata_combined = adata_combined[:, common_gene]
adata_combined.layers["spliced"] = adata_total.layers["spliced"]
adata_combined.layers["unspliced"] = adata_total.layers["unspliced"]
output_file = os.path.join(data_path, 'scvelo_merged.h5ad')
adata_combined.write(output_file)