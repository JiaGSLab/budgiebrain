#https://cole-trapnell-lab.github.io/monocle3/docs/introduction/
library(Seurat)
library(monocle3)
library(scales)
library(ggplot2)
library(dplyr)
budgerigar_dev_Tel_Glu <- readRDS('/data/work/Develop/12_budgerigar_dev_Tel_Glu_harmony.rds')
expression_matrix = budgerigar_dev_Tel_Glu@assays$RNA@counts
cell_metadata = data.frame(budgerigar_dev_Tel_Glu@meta.data)
gene_annotation = data.frame(gene_short_name = row.names(budgerigar_dev_Tel_Glu))
rownames(gene_annotation) =row.names(budgerigar_dev_Tel_Glu)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "Sample_ID")
## Step 3:cluster
cds <- reduce_dimension(cds,cores=5)
cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- learn_graph(cds)
## Step 4:original Embeddings
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(budgerigar_dev_Tel_Glu, reduction = "umap_harmony")
int.embed <- int.embed[rownames(cds.embed),]
## Step 5:calculate pseudotime
myselect <- function(cds,select.classify,my_select){ cell_ids <- which(colData(cds)[,select.classify] == my_select)  
                                                    closest_vertex <-    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex  
                                                    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])  
                                                    root_pr_nodes <-    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names    
                                                                                                                                  (which.max(table(closest_vertex[cell_ids,]))))]  
                                                    root_pr_nodes}
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'harmony_snn_res.0.9',my_select = "0"))
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, color_cells_by = "pseudotime",           
           show_trajectory_graph=F) 

