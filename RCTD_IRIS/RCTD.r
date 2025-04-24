library(Seurat)
library(dplyr)
downsampled_cells <- c()
annotations <- unique(budgerigar_5months_Tel@meta.data$subcelltype_annotation)
Idents(budgerigar_5months_Tel) <- budgerigar_5months_Tel$`subcelltype_annotation`
for (annotation in annotations) {
    cells_in_cluster <- WhichCells(budgerigar_5months_Tel, ident = annotation)
    if (length(cells_in_cluster) > 100) {
        set.seed(42) 
        sampled_cells <- sample(cells_in_cluster, 100)
    } else {
        sampled_cells <- cells_in_cluster
    }
    downsampled_cells <- c(downsampled_cells, sampled_cells)
}
saveRDS(seurat_downsampled, file = '/data/work/5month_Tel/5months_Tel_downsampled.rds')

# 使用下采样后的细胞创建新的Seurat对象
seurat_downsampled <- subset(budgerigar_5months_Tel, cells = downsampled_cells)
sc_input <- readRDS('/data/work/5month_Tel/5months_Tel_downsampled.rds')
sc_input <- subset(sc_input, subset = subclass_annotation != "Blood")
sc_input <- subset(sc_input, subset = subclass_annotation != "Mural cell")
sc_input <- subset(sc_input, subset = subclass_annotation != "Astrocyte")
sc_input <- subset(sc_input, subset = subclass_annotation != "OPC")
sc_input <- subset(sc_input, subset = subclass_annotation != "Microglia")
sc_input <- subset(sc_input, subset = subclass_annotation != "Astroependymal")
sc_input <- subset(sc_input, subset = subclass_annotation != "Oligodendrocyte")
sc_count <- sc_input@assays$RNA@counts
sc_meta <- sc_input@meta.data
dge_df <- as.data.frame(as.matrix(sc_count))
write.csv(dge_df, file = "/data/work/RCTD/Tel_scRNA/dge.csv", row.names = TRUE)
#subclasses
meta_data_df <- data.frame(
  barcode = rownames(sc_meta),
  cluster = sc_meta$subclass_annotation, 
  nUMI = sc_meta$nCount_RNA 
)

# csv_subclasses
write.csv(meta_data_df, file = "/data/work/RCTD/Tel_scRNA/meta_data_subclass.csv", row.names = FALSE)
#supertypes
meta_data_df <- data.frame(
  barcode = rownames(sc_meta),
  cluster = sc_meta$subcelltype_annotation, 
  nUMI = sc_meta$nCount_RNA 
)

# csv_supertypes
write.csv(meta_data_df, file = "/data/work/RCTD/Tel_scRNA/meta_data_subcelltype.csv", row.names = FALSE)
sp_input <- readRDS("/data/work/RCTD/5mm/5mm_raw!Vasc.rds")
meta_data <- sp_input@meta.data
spatial_location <-meta_data %>% 
    select(x,y)
bead_locations <- data.frame(
  barcodes = rownames(spatial_location), 
  xcoord = spatial_location[, "x"], 
  ycoord = spatial_location[, "y"]  
)
sp_count <- sp_input@assays$RNA@counts
mapped_dge_df <- as.data.frame(as.matrix(sp_count))
write.csv(mapped_dge_df, file = "/data/work/RCTD/5mm/MappedDGEForR.csv", row.names = TRUE)
write.csv(bead_locations, file = "/data/work/RCTD/5mm/BeadLocationsForR.csv", row.names = FALSE)



library(spacexr)
library(Matrix)
### Load in/preprocess your data, this might vary based on your file type
#refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
#counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
counts <- read.csv("/data/work/RCTD/Tel_scRNA/dge.csv")
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
#supertype
meta_data <- read.csv("/data/work/RCTD/Tel_scRNA/meta_data_subcelltype.csv") # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
table(reference@cell_types) #number of occurences for each cell type
saveRDS(reference, '/data/work/RCTD/5mm/SCRef.rds')
counts <- read.csv("/data/work/RCTD/5mm/MappedDGEForR.csv") # load in counts matrix
coords <- read.csv("/data/work/RCTD/5mm/BeadLocationsForR.csv")
coords$barcodes = paste0("X", coords$barcodes)
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL
head(coords)
rownames(counts) <- counts[,1]
counts[,1] <- NULL
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). # This list can be restricted if you want to crop the puck e.g. # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI')
#run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
norm_weights_matrix <- as.matrix(norm_weights)
# save
write.csv(norm_weights_matrix, file = "/data/work/RCTD/5mm/5mm_norm_weights.csv", row.names = TRUE)