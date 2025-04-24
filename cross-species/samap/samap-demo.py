from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
A=pd.read_csv('/data/work/samap/Tel/protein/maps/BGMM/BG_to_MM.txt',sep='\t',index_col=0,header=None)
B=pd.read_csv('/data/work/samap/Tel/protein/maps/BGMM/MM_to_BG.txt',sep='\t',index_col=0,header=None)
fn1 = '/data/work/samap/Tel/downsampled_GABA.h5ad'
fn2 = '/data/work/samap/cross-species/mousebrain/yao_2021_gaba_adata.h5ad'
filenames = {'BG':fn1,'MM':fn2}
sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)
sam1.adata.var_names=["BG_"+idx for idx in sam1.adata.var_names]
sam2.adata.var_names=["MM_"+idx for idx in sam2.adata.var_names]
sam1.adata.var_names
sam2.adata.var_names
import numpy as np
sam1.adata.var["weights"] = np.ones(sam1.adata.shape[1])
sam2.adata.var["weights"] = np.ones(sam2.adata.shape[1])


import scanpy as sc
import numpy as np
sam1 = SAM()
sam1.load_data(fn1) 
sc.pp.neighbors(sam1.adata)  
sam1.adata.uns["run_args"] = {"preprocessing": "StandardScaler", "weight_PCs": False}  
sam1.adata.var["weights"] = np.ones(sam1.adata.shape[1]) 

sam2 = SAM()
sam2.load_data(fn2)
sc.pp.neighbors(sam2.adata)  
sam2.adata.uns["run_args"] = {"preprocessing": "StandardScaler", "weight_PCs": False}  
sam2.adata.var["weights"] = np.ones(sam2.adata.shape[1])  

sams = {'BG': sam1, 'MM': sam2}

#create SAMap object
sm = SAMAP(
    sams,
    f_maps='/data/work/samap/Tel/protein/maps/',
)

sm.run(pairwise=True)
samap = sm.samap 
#定义keys字典：选取物种聚类标签。在这个例子中，'MM'和'TG'物种的聚类标签列名都是'leiden_clusters'。
keys = {'BG':'GABA_level_subcelltype_annotation','MM':'supertype_label'}

#计算映射得分并生成映射表：
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)#调用get_mapping_scores函数，传递sm对象、keys字典和参数n_top=0

MappingTable.to_csv('/data/work/samap/MappingTable_BG_MM_2021_GABA_supertypes.csv', index=True, header=True)

#定义keys字典：选取物种聚类标签。在这个例子中，'MM'和'TG'物种的聚类标签列名都是'leiden_clusters'。
keys = {'BG':'GABA_level_subclass_annotation','MM':'subclass_label'}

#计算映射得分并生成映射表：
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)#调用get_mapping_scores函数，传递sm对象、keys字典和参数n_top=0
#sm：包含多个物种的SAM对象。
#keys：包含每个物种的聚类标签。
#n_top=0：指定要返回的顶级映射得分的数量。n_top=0通常表示返回所有映射得分。
#######返回值：
#D：这是一个矩阵，包含了不同物种之间的映射得分。每个元素表示两个物种之间的映射得分。
#MappingTable：这是一个表格，包含了详细的映射信息，例如每个物种的细胞类型和对应的映射得分。

MappingTable.to_csv('/data/work/samap/MappingTable_BG_MM_2021_GABA_subclasses.csv', index=True, header=True)