import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse
# path
path = '/data/input/Files/RNA_velo/Tel_Dev/B01DA/4265-1-220613/' 
barcodes = pd.read_csv(path + 'barcodes.tsv.gz', header=None).to_numpy().flatten()
features = pd.read_csv(path + 'features.tsv.gz', header=None).to_numpy().flatten()

#  spliced and unspliced 
spliced = scipy.io.mmread(path + 'spliced.mtx.gz').T.tocsr()  
unspliced = scipy.io.mmread(path + 'unspliced.mtx.gz').T.tocsr()  

# create object
adata = sc.AnnData(X=spliced)
adata.layers['spliced'] = spliced
adata.layers['unspliced'] = unspliced

# set name
adata.obs_names = barcodes
adata.var_names = features

# align
assert adata.shape[0] == len(adata.obs_names), "cell numbers wrong"
assert adata.shape[1] == len(adata.var_names), "gene counts wrong"

print(adata)
sc.write('/data/work/Tel_dev/scVelo/data/B01DA_1',adata)
