library(sceasy)
library(reticulate)
h5ad_file <-"/data/work/figures/annotation/5mm.h5ad"
sceasy::convertFormat(file.path(h5ad_file), from="anndata", to="seurat",
                       outFile=file.path('/data/work/IRIS/5mm_BS16_iris.rds'))




library(Seurat)
library(dplyr)
#downsampled snRNAseq
sc_input <- readRDS('/data/work/5month_Tel/5months_Tel_downsampled.rds')
sc_count <- sc_input@assays$RNA@counts
sc_meta <- sc_input@meta.data
ct.varname <- "subcelltype_annotation"
sample.varname <- "Sample_ID"
#smoothed ST
BS16_5mm <-readRDS("/data/work/IRIS/5mm_BS16_iris.rds")
sp_input <- subset(BS16_5mm, subset = simplified_annotated_spatial_leiden != "Vasc")
meta_data <- sp_input@meta.data
spatial_location_list <- list(meta_data %>% 
    select(x,y))
#raw data
sp_input2 <- readRDS("/data/work/converted_rds/C02946C3.tissue.rds")
cell_names <- rownames(spatial_location_list[[1]])
sp_input2_subset <- subset(sp_input2, cells = cell_names)
print(sp_input2_subset)
spatial_countMat_list <- list(sp_input2_subset@assays$RNA@counts)

#IRIS predict cell distribution
library(IRIS)
IRIS_object <- createIRISObject(
spatial_countMat_list = spatial_countMat_list,
spatial_location_list = spatial_location_list,
sc_count = sc_count,
sc_meta = sc_meta,
ct.varname = ct.varname,
sample.varname = sample.varname,
minCountGene = 100,
minCountSpot =5) 

numCluster =6 #determined by brain regions
IRIS_object <- IRIS_spatial(IRIS_object,numCluster = numCluster)

colors = c("#ebe5c2", "#D57358", "#023047", "#F7CA71", "#1697a6", "#8bc6cc", "#C9DEC3")
#### extract the domain labels detected by IRIS
IRIS_domain =  IRIS_object@spatialDomain[,c("Slice","spotName","IRIS_domain")]
#### extract the spatial location information
spatial_location = IRIS_object@spatialDomain[,c("x","y")]
spatial_location$x = IRIS_object@spatialDomain$y
spatial_location$y = -IRIS_object@spatialDomain$x 
p1 = IRIS.visualize.domain(IRIS_domain,spatial_location,colors = colors,numCols = 1)
print(p1)
#### select the cell type that we are interested
ct.visualize = c("DACH2_NDST4_Ex","DACH2_MEIS2_Ex","DACH2_LYPD1_Ex","DACH2_BDNF_Ex","DACH2_ZNF381_Ex","DACH2_CHD7_Ex")
#### visualize the spatial distribution of the cell type proportion
library(viridis)
#### set the color values
colors = magma(256)
#### extract the cell type proportion matrix for the example slice
IRIS_prop = IRIS_object@IRIS_Prop
IRIS_prop = IRIS_prop[IRIS_prop$Slice == "Slice1",]
#### extract the spatial information matrix for the example slice
IRIS_prop_location = IRIS_prop[,c("Slice","spotName","x","y")]
#### This is only for visualization purpose
IRIS_prop_location$x = IRIS_prop$y
IRIS_prop_location$y = -IRIS_prop$x
p3 <- IRIS.visualize.eachProp(
	proportion = IRIS_prop,        
	spatial_location = IRIS_prop_location, 
	ct.visualize = ct.visualize,                 #### selected cell types to visualize
	colors = colors,                             #### if not provide, we will use the default colors
	numCols = 2)                                 #### number of columns in the figure panel

ggsave("/data/work/IRIS/Tel/5mm/Tel_celltype2.pdf", width = 6, height = 7)